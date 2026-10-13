"""
sp_tags pipeline: long reads -> split molecules -> spatial_id / molecule_id -> tagged FASTQ.

Step 1  Detect tags in every read (in parallel); concatenated reads are split into
        molecules. Writes <prefix>.tags.tsv and a temporary per-molecule FASTQ.
Step 2  Correct molecule_id (UMI) errors within each spatial_id and write the final
        <prefix>.fastq, ready for minimap2 and `main.py quantify --st_mode`.

Quantification unit (what is written as cell barcode CB):
    default          CB = spatial_id (CID),  UB = molecule_id (MID)
    --barcode2spot   CB = spot / cell the CID belongs to,  UB = CID + MID
                     (the CID keeps molecules from different CIDs of one spot apart
                     during UMI deduplication); CIDs absent from the mapping are dropped.

Read headers of the final FASTQ:
    @<segment_id>_<CB>_<UB>\tCB:Z:<CB>\tCR:Z:<raw CID>\tUB:Z:<UB>\tUR:Z:<raw MID>[\tXC:Z:<CID>]
(XC only with --barcode2spot), so both barcode modes of miniTEQuant work:
    --barcode_in_readname       (default)  parses <CB>_<UB> from the read name
    --no_barcode_in_readname               reads CB / UB tags (align with `minimap2 -y`)
"""

import argparse
import multiprocessing as mp
import os
import shutil
import time
from collections import defaultdict, deque

from .align import BACKEND
from .base import TagResult, TagStats
from .kmer_index import BarcodeIndex
from .seq_utils import detect_format, iter_reads, reverse_complement
from .spot_map import load_spot_map, parse_spec
from .stereo_seq import StereoSeqDetector
from .umi import cluster_umis
from .whitelist import load_whitelist

PLATFORMS = {
    'stereo': StereoSeqDetector,
}

# Read technology -> defaults. Adapter detection is the same for both (as in spl-IsoQuant);
# nanopore reads need a shorter whitelist k-mer to retrieve CIDs carrying several errors.
DATA_TYPES = {'pacbio': 'pacbio', 'nanopore': 'nanopore', 'ont': 'nanopore'}
KMER_SIZE = {'pacbio': 14, 'nanopore': 12}

# Detector shared with worker processes (inherited on fork, set by _init_worker on spawn).
_worker_detector = None
_worker_options = None


def log(msg):
    print('[sp_tags] %s' % msg, flush=True)


# ---------------------------------------------------------------------------
# Step 1: tag detection
# ---------------------------------------------------------------------------

class _Options:
    def __init__(self, split, require_good_umi, fastq_output):
        self.split = split
        self.require_good_umi = require_good_umi
        self.fastq_output = fastq_output


def _init_worker(detector, options):
    global _worker_detector, _worker_options
    _worker_detector = detector
    _worker_options = options


def _format_record(fastq_output, name, comment, seq, qual):
    if fastq_output:
        if qual is None:
            qual = '!' * len(seq)
        return '@%s\t%s\n%s\n+\n%s\n' % (name, comment, seq, qual)
    return '>%s\t%s\n%s\n' % (name, comment, seq)


def _process_read(detector, options, read_id, seq, qual, stats, tsv_lines, records):
    results = detector.detect(read_id, seq)
    stats.add('Input reads')

    # In a concatenated read only molecules closed by a TSO are complete.
    require_tso = options.split and len(results) > 1
    rc_seq, rc_qual = None, None
    n_valid = 0
    for r in results:
        stats.add_result(r)
        if not r.is_valid():
            tsv_lines.append(r.to_tsv())
            continue
        n_valid += 1

        if r.strand == '-':
            if rc_seq is None:
                rc_seq = reverse_complement(seq)
                rc_qual = qual[::-1] if qual is not None else None
            oriented_seq, oriented_qual = rc_seq, rc_qual
        else:
            oriented_seq, oriented_qual = seq, qual

        if options.split:
            start, end = r.segment_start(), r.segment_end(len(seq))
            r.read_id = '%s_%d_%d_%s' % (read_id, start, end, r.strand)
        else:
            start, end = 0, len(seq)
        tsv_lines.append(r.to_tsv())

        if require_tso and r.tso5 == -1:
            stats.add('Molecules skipped: no TSO in concatenated read')
            continue
        if not r.has_umi():
            stats.add('Molecules skipped: no UMI')
            continue
        if options.require_good_umi and not r.umi_good:
            stats.add('Molecules skipped: UMI length out of range')
            continue

        comment = '%s\t%s\t%s\t%d' % (r.spatial_id, r.raw_spatial_id, r.molecule_id, r.whitelist_index)
        segment_qual = oriented_qual[start:end] if oriented_qual is not None else None
        records.append(_format_record(options.fastq_output, r.read_id, comment,
                                      oriented_seq[start:end], segment_qual))
        stats.add('Molecules written')

    if options.split:
        if n_valid == 0:
            stats.add('Reads with 0 molecules')
        elif n_valid == 1:
            stats.add('Reads with 1 molecule')
        elif n_valid == 2:
            stats.add('Reads with 2 molecules')
        else:
            stats.add('Reads with 3+ molecules')


def _process_chunk(chunk):
    stats = TagStats()
    tsv_lines, records = [], []
    for read_id, seq, qual in chunk:
        if seq:
            _process_read(_worker_detector, _worker_options, read_id, seq, qual, stats, tsv_lines, records)
    tsv_text = '\n'.join(tsv_lines) + '\n' if tsv_lines else ''
    return tsv_text, ''.join(records), dict(stats.counts), len(chunk)


def _read_chunks(input_files, chunk_size):
    chunk = []
    for path in input_files:
        log('Reading %s' % path)
        for read in iter_reads(path):
            chunk.append(read)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    if chunk:
        yield chunk


def detect_tags(input_files, detector, options, tsv_path, tmp_records_path, threads, chunk_size):
    stats = TagStats()
    processed = [0]

    with open(tsv_path, 'w') as tsv_out, open(tmp_records_path, 'w') as rec_out:
        tsv_out.write(TagResult.header(detector.SPATIAL_LABEL, detector.MOLECULE_LABEL) + '\n')

        def write(result):
            tsv_text, records_text, counts, n_reads = result
            tsv_out.write(tsv_text)
            rec_out.write(records_text)
            stats.merge(counts)
            processed[0] += n_reads
            log('Processed %d reads' % processed[0])

        chunks = _read_chunks(input_files, chunk_size)
        if threads <= 1:
            _init_worker(detector, options)
            for chunk in chunks:
                write(_process_chunk(chunk))
            return stats

        if 'fork' in mp.get_all_start_methods():
            # workers inherit the (possibly very large) whitelist index copy-on-write
            _init_worker(detector, options)
            pool = mp.get_context('fork').Pool(threads)
        else:
            pool = mp.get_context('spawn').Pool(threads, initializer=_init_worker, initargs=(detector, options))

        try:
            pending = deque()
            for chunk in chunks:
                pending.append(pool.apply_async(_process_chunk, (chunk,)))
                # keep a bounded number of chunks in flight, and write them in input order
                while len(pending) >= 2 * threads:
                    write(pending.popleft().get())
            while pending:
                write(pending.popleft().get())
        finally:
            pool.close()
            pool.join()
    return stats


# ---------------------------------------------------------------------------
# Step 2: UMI correction and final output
# ---------------------------------------------------------------------------

def _iter_records(path, fastq_output):
    lines_per_record = 4 if fastq_output else 2
    with open(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            body = [handle.readline() for _ in range(lines_per_record - 1)]
            fields = header[1:].rstrip('\n').split('\t')
            yield fields, body


def count_umis(tmp_records_path, fastq_output):
    """Reads per (spatial_id, molecule_id)."""
    umi_counts = defaultdict(lambda: defaultdict(int))
    for fields, _ in _iter_records(tmp_records_path, fastq_output):
        umi_counts[fields[1]][fields[3]] += 1
    return umi_counts


def correct_and_write(tmp_records_path, out_path, fastq_output, umi_counts, umi_ed, umi_length, name_sep,
                      stats, spot_map=None):
    """
    Correct UMIs within each spatial_id and write the tagged reads.

    spot_map : {spatial_id: spot} from --barcode2spot, or None to use the spatial_id as CB.
    Returns ({spatial_id: whitelist index}, {CB: reads written}).
    """

    corrections = {}
    n_umis, n_molecules = 0, 0
    for spatial_id, counts in umi_counts.items():
        mapping = cluster_umis(counts, umi_ed, umi_length)
        n_umis += len(mapping)
        n_molecules += len(set(mapping.values()))
        for umi, representative in mapping.items():
            if umi != representative:
                corrections[(spatial_id, umi)] = representative
    stats.add('Spatial IDs with molecules', len(umi_counts))
    stats.add('Distinct raw UMIs (per spatial ID)', n_umis)
    stats.add('Distinct corrected UMIs (per spatial ID)', n_molecules)
    del umi_counts

    used_whitelist_index = {}
    cb_reads = defaultdict(int)
    marker = '@' if fastq_output else '>'
    with open(out_path, 'w') as out:
        for fields, body in _iter_records(tmp_records_path, fastq_output):
            segment_id, spatial_id, raw_spatial_id, raw_umi, wl_index = fields
            umi = corrections.get((spatial_id, raw_umi), raw_umi)
            if spot_map is None:
                cb, ub, extra = spatial_id, umi, ''
            else:
                cb = spot_map.get(spatial_id)
                if cb is None:
                    stats.add('Molecules dropped: spatial ID not in barcode2spot')
                    continue
                # the CID keeps equal MIDs of different CIDs in one spot apart
                ub, extra = spatial_id + umi, '\tXC:Z:' + spatial_id
            used_whitelist_index[spatial_id] = int(wl_index)
            cb_reads[cb] += 1
            out.write('%s%s%s%s%s%s\tCB:Z:%s\tCR:Z:%s\tUB:Z:%s\tUR:Z:%s%s\n' % (
                marker, segment_id, name_sep, cb, name_sep, ub,
                cb, raw_spatial_id, ub, raw_umi, extra))
            out.writelines(body)
    if spot_map is not None:
        stats.add('Molecules assigned to spots', sum(cb_reads.values()))
        stats.add('Spatial IDs assigned to spots', len(used_whitelist_index))
        stats.add('Spots with molecules', len(cb_reads))
    return used_whitelist_index, cb_reads


def write_spatial_positions(path, used_whitelist_index, whitelist):
    """Write detected CIDs with coordinates in tissue_positions.csv format (for --tissue_positions)."""
    with open(path, 'w') as out:
        out.write('barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n')
        for spatial_id in sorted(used_whitelist_index):
            i = used_whitelist_index[spatial_id]
            if i < 0:
                continue
            x, y = int(whitelist.x[i]), int(whitelist.y[i])
            out.write('%s,1,%d,%d,%d,%d\n' % (spatial_id, y, x, y, x))


def write_spot_summary(summary_path, positions_path, used_whitelist_index, spot_map, cb_reads, whitelist):
    """
    Per-spot summary (spot, spatial IDs, reads, centroid) and, if the whitelist has
    coordinates, spot centroids in tissue_positions.csv format.
    """
    members = defaultdict(list)
    for spatial_id, i in used_whitelist_index.items():
        members[spot_map[spatial_id]].append(i)
    with_coords = whitelist is not None and whitelist.has_coordinates

    def centroid(indices):
        indices = [i for i in indices if i >= 0]
        if not with_coords or not indices:
            return None
        return (int(round(sum(int(whitelist.x[i]) for i in indices) / len(indices))),
                int(round(sum(int(whitelist.y[i]) for i in indices) / len(indices))))

    centroids = {spot: centroid(indices) for spot, indices in members.items()}
    with open(summary_path, 'w') as out:
        out.write('spot\tn_spatial_ids\tn_reads\tx_center\ty_center\n')
        for spot in sorted(members):
            c = centroids[spot]
            out.write('%s\t%d\t%d\t%s\t%s\n' % (spot, len(members[spot]), cb_reads[spot],
                                                c[0] if c else 'NA', c[1] if c else 'NA'))
    if not with_coords:
        return
    with open(positions_path, 'w') as out:
        out.write('barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n')
        for spot in sorted(members):
            c = centroids[spot]
            if c is not None:
                out.write('%s,1,%d,%d,%d,%d\n' % (spot, c[1], c[0], c[1], c[0]))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run(args):
    start_time = time.time()
    os.makedirs(args.output_dir, exist_ok=True)
    prefix = os.path.join(args.output_dir, args.prefix)
    tmp_dir = os.path.join(args.output_dir, 'tmp_%s' % args.prefix)
    os.makedirs(tmp_dir, exist_ok=True)

    if not args.whitelist:
        log('WARNING: no whitelist given; raw (uncorrected) spatial IDs will be reported.')
    if BACKEND == 'python':
        log('WARNING: neither ssw-py nor numba is installed; using slow pure-Python alignment.')

    if args.barcode2spot:
        parse_spec(args.barcode2spot)    # fail early on a missing file
        log('Quantification unit: spot (spatial ID -> spot mapping from %s)' % args.barcode2spot)
    else:
        log('Quantification unit: spatial ID (no --barcode2spot given)')

    input_formats = {detect_format(p) for p in args.input}
    fastq_output = 'fasta' not in input_formats
    extension = 'fastq' if fastq_output else 'fasta'

    data_type = DATA_TYPES[args.data_type]
    if args.kmer_size is None:
        args.kmer_size = KMER_SIZE[data_type]
    log('Data type: %s' % data_type)

    whitelist, barcode_index = None, None
    if args.whitelist:
        log('Loading whitelist %s' % args.whitelist)
        whitelist = load_whitelist(args.whitelist)
        log('Loaded %d barcodes of length %d%s' % (len(whitelist), whitelist.seq_len,
                                                    ' with coordinates' if whitelist.has_coordinates else ''))
        barcode_index = BarcodeIndex(whitelist.barcodes, whitelist.seq_len, kmer_size=args.kmer_size)
        log('Built %s k-mer index (k=%d)' % ('dense' if barcode_index.dense else 'sparse', args.kmer_size))

    detector_class = PLATFORMS[args.platform]
    detector = detector_class(barcode_index, split=not args.no_split, min_score=args.min_bc_score,
                              primer=args.primer, linker=args.linker, tso=args.tso)
    log('Adapters: primer=%s linker=%s TSO=%s' % (detector.primer, detector.linker, detector.tso))
    options = _Options(split=not args.no_split, require_good_umi=args.require_good_umi,
                       fastq_output=fastq_output)

    tsv_path = prefix + '.tags.tsv'
    tmp_records_path = os.path.join(tmp_dir, 'molecules.' + extension)
    log('Step 1/2: detecting %s / %s (alignment backend: %s, threads: %d)' % (
        detector.SPATIAL_LABEL, detector.MOLECULE_LABEL, BACKEND, args.threads))
    stats = detect_tags(args.input, detector, options, tsv_path, tmp_records_path,
                        args.threads, args.chunk_size)

    log('Step 2/2: correcting %s errors (max edit distance %d) and writing output' % (
        detector.MOLECULE_LABEL, args.umi_ed))
    out_reads_path = '%s.%s' % (prefix, extension)
    umi_counts = count_umis(tmp_records_path, fastq_output)

    spot_map = None
    if args.barcode2spot:
        label = detector.SPATIAL_LABEL
        spot_map, info = load_spot_map(args.barcode2spot, args.name_sep, keep_ids=set(umi_counts))
        log('barcode2spot: %d of %d detected %ss mapped to %d spots (columns %d:%d)' % (
            len(spot_map), len(umi_counts), label, len(set(spot_map.values())), *info['columns']))
        if info['skipped']:
            log('barcode2spot: skipped %d header / malformed lines' % info['skipped'])
        if info['conflicts']:
            log('WARNING: barcode2spot: %d %ss map to several spots; the first one is used'
                % (info['conflicts'], label))
        if info['renamed']:
            old, new = next(iter(info['renamed'].items()))
            log("WARNING: barcode2spot: %d spot names changed to be valid cell barcodes, e.g. '%s' -> '%s'"
                % (len(info['renamed']), old, new))
        if not spot_map:
            log('WARNING: barcode2spot: no detected %s is in the mapping file; check that column %d '
                'holds %s sequences' % (label, info['columns'][0], label))

    used_whitelist_index, cb_reads = correct_and_write(
        tmp_records_path, out_reads_path, fastq_output, umi_counts, args.umi_ed, detector.UMI_LEN,
        args.name_sep, stats, spot_map)

    if spot_map is not None:
        summary_path = prefix + '.spots.tsv'
        positions_path = prefix + '.spot_positions.csv'
        write_spot_summary(summary_path, positions_path, used_whitelist_index, spot_map, cb_reads, whitelist)
        log('Spot summary: %s' % summary_path)
        if whitelist is not None and whitelist.has_coordinates:
            log('Spot positions (centroids): %s' % positions_path)
    elif whitelist is not None and whitelist.has_coordinates:
        positions_path = prefix + '.spatial_positions.csv'
        write_spatial_positions(positions_path, used_whitelist_index, whitelist)
        log('Spatial positions: %s' % positions_path)

    stats_path = prefix + '.stats.tsv'
    with open(stats_path, 'w') as out:
        out.write('\n'.join(stats.lines()) + '\n')
    for line in stats.lines():
        log('  ' + line.replace('\t', ': '))

    if not args.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    log('Tag table   : %s' % tsv_path)
    log('Tagged reads: %s' % out_reads_path)
    log('Statistics  : %s' % stats_path)
    log('Finished in %.1f min' % ((time.time() - start_time) / 60))


def build_parser():
    parser = argparse.ArgumentParser(
        description='Detect spatial IDs (barcode / CID) and molecule IDs (UMI / MID) in long reads, '
                    'split concatenated reads into molecules and write tagged reads for miniTEQuant.')
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='Long-read files: FASTQ / FASTA (optionally .gz) or unaligned BAM')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('--prefix', default='sp_tags', help='Output file prefix (default: sp_tags)')
    parser.add_argument('--platform', choices=sorted(PLATFORMS), default='stereo',
                        help='Spatial platform (default: stereo)')
    parser.add_argument('-w', '--whitelist', default=None,
                        help='Spatial barcode whitelist: one barcode per line, optionally followed by '
                             'x and y coordinates (e.g. Stereo-seq "CID x y"); plain text or .gz')
    parser.add_argument('--barcode2spot', default=None, metavar='FILE[:ID_COL:SPOT_COL]',
                        help='Mapping of spatial IDs (CIDs) to spots / cells, e.g. a cellbin assignment: '
                             'spatial ID in column 0 and spot in column 1 by default, or 0-based columns '
                             'given as FILE:ID_COL:SPOT_COL. Reads are then tagged with CB = spot and '
                             'UB = CID+MID so quantification is per spot; CIDs absent from the file are '
                             'dropped. Without it, quantification is per CID.')
    parser.add_argument('--no_split', action='store_true', default=False,
                        help='Do not split concatenated reads; report one molecule per read')
    parser.add_argument('--primer', default=None,
                        help='Primer upstream of the spatial ID, 5\'->3\' in capture-oligo-first orientation '
                             '(default for stereo: %s)' % StereoSeqDetector.PC1_PRIMER)
    parser.add_argument('--linker', default=None,
                        help='Linker between spatial ID and UMI (default for stereo: %s)' % StereoSeqDetector.LINKER)
    parser.add_argument('--tso', default=None,
                        help='TSO marking the end of a molecule (default for stereo: %s)' % StereoSeqDetector.TSO5)
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of processes (default: 1)')
    parser.add_argument('--umi_ed', type=int, default=1,
                        help='Max edit distance to merge UMIs within a spatial ID; 0 disables correction '
                             '(default: 1)')
    parser.add_argument('--require_good_umi', action='store_true', default=False,
                        help='Only write molecules whose UMI length is within the expected range')
    parser.add_argument('--min_bc_score', type=int, default=None,
                        help='Min alignment score to accept a whitelist barcode (default: barcode length - 3)')
    parser.add_argument('--data_type', choices=sorted(DATA_TYPES), default='pacbio',
                        help='Long-read technology; sets defaults for higher-error reads. nanopore (alias ont): '
                             '--kmer_size 12, so that whitelist CIDs are still found when sequencing errors '
                             'destroy most 14-mers (default: pacbio)')
    parser.add_argument('--kmer_size', type=int, default=None,
                        help='K-mer size of the whitelist index (default: 14 for pacbio, 12 for nanopore)')
    parser.add_argument('--name_sep', default='_',
                        help="Separator between read name, barcode and UMI in output read names; must match "
                             "--barcode_separator of main.py quantify (default: '_')")
    parser.add_argument('--chunk_size', type=int, default=10000,
                        help='Reads per parallel work unit (default: 10000)')
    parser.add_argument('--keep_tmp', action='store_true', default=False, help='Keep temporary files')
    return parser


def main(argv=None):
    run(build_parser().parse_args(argv))
