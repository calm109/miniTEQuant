"""
Post-alignment UMI deduplication by gene (spl-IsoQuant style).

sp_tags corrects molecule_id (MID / UMI) errors before alignment, within each
spatial_id only, and therefore conservatively (--umi_ed 1). Long reads, ONT in
particular, carry more UMI errors than that. As in spl-IsoQuant
(barcode_calling/umi_filtering.py, UMIFilter), this step groups aligned reads by
(spatial_id, gene) and clusters their UMIs with a larger edit distance (Stereo-seq
default: 4): within one gene of one ~0.5 um spot, two genuine molecules rarely share
a UMI that similar. One representative read is kept per molecule:

    1. read assigned to a single gene  >  ambiguous assignment
    2. more exons (introns + 1)
    3. longer aligned length
    4. higher MAPQ

Reads overlapping no gene are kept unchanged (they are not quantified by miniTEQuant).

Input  : SAM / BAM of reads from sp_tags aligned with minimap2 (`-y` keeps the tags;
         without it, CB / UB are parsed from the read name). BAM needs samtools.
Output : SAM with one primary alignment per molecule (header kept), plus statistics.
"""

import argparse
import bisect
import gzip
import os
import re
import subprocess
import sys
import time
from collections import defaultdict

from .umi import cluster_umis

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')
_GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
_SKIP_FLAGS = 0x4 | 0x100 | 0x800    # unmapped, secondary, supplementary


def log(msg):
    print('[umi_dedup] %s' % msg, flush=True)


# ---------------------------------------------------------------------------
# Gene annotation
# ---------------------------------------------------------------------------

class GeneIndex:
    """Exons per chromosome, sorted by start, for overlap queries."""

    def __init__(self, gtf_path):
        gene_exons = defaultdict(list)
        opener = gzip.open if gtf_path.endswith('.gz') else open
        with opener(gtf_path, 'rt') as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                fields = line.split('\t', 8)
                if len(fields) < 9 or fields[2] != 'exon':
                    continue
                match = _GENE_ID_RE.search(fields[8])
                if match:
                    gene_exons[(fields[0], fields[6], match.group(1))].append((int(fields[3]) - 1, int(fields[4])))

        # merge the exons of all transcripts of a gene, so overlaps are not counted twice
        exons = defaultdict(list)
        for (chrom, strand, gene), intervals in gene_exons.items():
            intervals.sort()
            cur_start, cur_end = intervals[0]
            for start, end in intervals[1:]:
                if start <= cur_end:
                    cur_end = max(cur_end, end)
                else:
                    exons[chrom].append((cur_start, cur_end, strand, gene))
                    cur_start, cur_end = start, end
            exons[chrom].append((cur_start, cur_end, strand, gene))

        self.starts, self.exons, self.max_len = {}, {}, {}
        for chrom, items in exons.items():
            items.sort()
            self.exons[chrom] = items
            self.starts[chrom] = [e[0] for e in items]
            self.max_len[chrom] = max(e[1] - e[0] for e in items)
        self.n_genes = len({gene for _, _, gene in gene_exons})

    def overlaps(self, chrom, blocks):
        """{(gene_id, strand): overlapping bases} for 0-based half-open aligned blocks."""
        starts = self.starts.get(chrom)
        if starts is None:
            return {}
        exons, max_len = self.exons[chrom], self.max_len[chrom]
        result = defaultdict(int)
        for b_start, b_end in blocks:
            i = bisect.bisect_left(starts, b_start - max_len)
            j = bisect.bisect_left(starts, b_end)
            for e_start, e_end, strand, gene in exons[i:j]:
                overlap = min(b_end, e_end) - max(b_start, e_start)
                if overlap > 0:
                    result[(gene, strand)] += overlap
        return result


def aligned_blocks(pos, cigar):
    """Aligned reference blocks (0-based, half-open), number of exons, aligned length."""
    blocks = []
    ref = pos
    block_start = pos
    aligned = 0
    n_exons = 1
    for length, op in _CIGAR_RE.findall(cigar):
        length = int(length)
        if op in 'M=X':
            ref += length
            aligned += length
        elif op == 'D':
            ref += length
        elif op == 'N':
            if ref > block_start:
                blocks.append((block_start, ref))
            ref += length
            block_start = ref
            n_exons += 1
    if ref > block_start:
        blocks.append((block_start, ref))
    return blocks, n_exons, aligned


def assign_gene(gene_index, chrom, pos, cigar, reverse):
    """
    Return (gene_id or None, unambiguous, n_exons, aligned_length).

    Reads from sp_tags start with the capture oligo, i.e. they are antisense to the
    mRNA, so a forward alignment points to a gene on the '-' strand. Genes on the
    expected strand are preferred; the other strand is used only if nothing overlaps.
    """
    blocks, n_exons, aligned = aligned_blocks(pos, cigar)
    overlaps = gene_index.overlaps(chrom, blocks)
    if not overlaps:
        return None, False, n_exons, aligned
    expected = '+' if reverse else '-'
    candidates = {g: o for (g, s), o in overlaps.items() if s == expected}
    if not candidates:
        candidates = {g: o for (g, s), o in overlaps.items()}
    best = max(candidates.values())
    best_genes = sorted(g for g, o in candidates.items() if o == best)
    return best_genes[0], len(candidates) == 1, n_exons, aligned


# ---------------------------------------------------------------------------
# SAM input
# ---------------------------------------------------------------------------

def open_alignments(path, threads=1):
    """Yield SAM text lines from a SAM(.gz) or BAM file."""
    if path.endswith('.bam'):
        proc = subprocess.Popen(['samtools', 'view', '-h', '-@', str(max(1, threads - 1)), path],
                                stdout=subprocess.PIPE, universal_newlines=True, bufsize=1 << 20)
        try:
            for line in proc.stdout:
                yield line
        finally:
            proc.stdout.close()
            if proc.wait() != 0:
                raise RuntimeError('samtools view failed on %s' % path)
    else:
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line in handle:
                yield line


def read_tags(fields, name_sep):
    """(unit, umi) of a SAM record: XC (or CB) and the MID part of UB."""
    tags = {}
    for field in fields[11:]:
        if field[2:5] == ':Z:' and field[:2] in ('CB', 'UB', 'XC'):
            tags[field[:2]] = field[5:].rstrip('\n')
    if 'CB' not in tags or 'UB' not in tags:
        parts = fields[0].rsplit(name_sep, 2)
        if len(parts) != 3:
            return None, None
        tags.setdefault('CB', parts[1])
        tags.setdefault('UB', parts[2])
    cid = tags.get('XC')
    if cid and tags['UB'].startswith(cid):
        # --barcode2spot output: CB = spot, UB = CID + MID -> dedup within the CID
        return cid, tags['UB'][len(cid):]
    return tags['CB'], tags['UB']


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------

def collect(path, gene_index, name_sep, threads, stats):
    """
    First pass: for each primary alignment record (by ordinal), its gene group and
    selection key. Returns {(unit, gene): [(umi, rank_key, ordinal), ...]}.
    """
    groups = defaultdict(list)
    ordinal = -1
    for line in open_alignments(path, threads):
        if line.startswith('@'):
            continue
        ordinal += 1
        fields = line.split('\t')
        flag = int(fields[1])
        if flag & _SKIP_FLAGS:
            stats['Records skipped (unmapped / secondary / supplementary)'] += 1
            continue
        stats['Primary alignments'] += 1
        unit, umi = read_tags(fields, name_sep)
        if unit is None:
            stats['Primary alignments without barcode / UMI (kept)'] += 1
            continue
        gene, unambiguous, n_exons, aligned = assign_gene(
            gene_index, fields[2], int(fields[3]) - 1, fields[5], bool(flag & 0x10))
        if gene is None:
            stats['Primary alignments outside genes (kept)'] += 1
            continue
        if not unambiguous:
            stats['Primary alignments with ambiguous gene'] += 1
        rank = (unambiguous, n_exons, aligned, int(fields[4]))
        groups[(unit, gene)].append((umi, rank, ordinal))
    return groups


def select_representatives(groups, edit_distance, umi_length, stats):
    """Ordinals of the alignments to drop (duplicates of a kept molecule)."""
    drop = set()
    for reads in groups.values():
        counts = defaultdict(int)
        for umi, _, _ in reads:
            counts[umi] += 1
        mapping = cluster_umis(counts, edit_distance, umi_length)
        best = {}
        for umi, rank, ordinal in reads:
            molecule = mapping[umi]
            current = best.get(molecule)
            if current is None or rank > current[0]:
                if current is not None:
                    drop.add(current[1])
                best[molecule] = (rank, ordinal)
            else:
                drop.add(ordinal)
        stats['Molecules (UMI clusters) in genes'] += len(best)
        stats['Reads in genes'] += len(reads)
    stats['Gene / spatial ID groups'] = len(groups)
    stats['Duplicate reads removed'] = len(drop)
    return drop


def write_output(path, out_path, drop, threads, stats):
    ordinal = -1
    with open(out_path, 'w') as out:
        for line in open_alignments(path, threads):
            if line.startswith('@'):
                out.write(line)
                continue
            ordinal += 1
            if ordinal in drop or int(line.split('\t', 2)[1]) & _SKIP_FLAGS:
                continue
            out.write(line)
            stats['Alignments written'] += 1


def run(args):
    start_time = time.time()
    stats = defaultdict(int)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    log('Loading genes from %s' % args.gtf)
    gene_index = GeneIndex(args.gtf)
    log('Loaded %d genes on %d chromosomes' % (gene_index.n_genes, len(gene_index.exons)))

    log('Pass 1/2: assigning reads to genes (%s)' % args.input)
    groups = collect(args.input, gene_index, args.name_sep, args.threads, stats)
    log('Clustering UMIs per spatial ID and gene (max edit distance %d)' % args.umi_ed)
    drop = select_representatives(groups, args.umi_ed, args.umi_length, stats)
    del groups

    log('Pass 2/2: writing %s' % args.output)
    write_output(args.input, args.output, drop, args.threads, stats)

    reads, molecules = stats['Reads in genes'], stats['Molecules (UMI clusters) in genes']
    if molecules:
        stats['Mean reads per molecule (x100)'] = int(round(100.0 * reads / molecules))
    lines = ['%s\t%d' % (k, stats[k]) for k in sorted(stats)]
    stats_path = args.stats or os.path.splitext(args.output)[0] + '.umi_dedup.stats.tsv'
    with open(stats_path, 'w') as out:
        out.write('\n'.join(lines) + '\n')
    for line in lines:
        log('  ' + line.replace('\t', ': '))
    log('Statistics: %s' % stats_path)
    log('Finished in %.1f min' % ((time.time() - start_time) / 60))


def build_parser():
    parser = argparse.ArgumentParser(
        description='Post-alignment UMI deduplication per spatial ID and gene for reads produced by '
                    'sp_tags (spl-IsoQuant style); keeps one alignment per molecule.')
    parser.add_argument('-i', '--input', required=True,
                        help='Alignments of sp_tags reads: SAM, SAM.gz or BAM (BAM needs samtools in PATH)')
    parser.add_argument('-g', '--gtf', required=True, help='Gene annotation (GTF, optionally .gz)')
    parser.add_argument('-o', '--output', required=True, help='Output SAM with one alignment per molecule')
    parser.add_argument('--umi_ed', type=int, default=4,
                        help='Max edit distance to merge UMIs within one spatial ID and gene '
                             '(spl-IsoQuant Stereo-seq default: 4)')
    parser.add_argument('--umi_length', type=int, default=10,
                        help='Expected UMI length, used to break ties between UMIs (default: 10, Stereo-seq MID)')
    parser.add_argument('--name_sep', default='_',
                        help="Read-name separator used by sp_tags, for alignments without CB/UB tags (default: '_')")
    parser.add_argument('--stats', default=None, help='Statistics file (default: <output>.umi_dedup.stats.tsv)')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Threads for samtools when reading BAM')
    return parser


def main(argv=None):
    run(build_parser().parse_args(argv))


if __name__ == '__main__':
    sys.exit(main())
