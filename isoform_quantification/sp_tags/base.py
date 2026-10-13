"""
Platform-independent result containers for spatial tag detection.

Every platform maps its own identifiers onto two unified fields:
    spatial_id  : position identifier (Stereo-seq CID, Visium HD barcode, ...)
    molecule_id : molecule identifier (Stereo-seq MID, 10x UMI, ...)

Ported from spl-IsoQuant barcode_calling/callers/base.py
(LinkerBarcodeDetectionResult / TSOBarcodeDetectionResult / SplittingBarcodeDetectionResult).
"""

from collections import defaultdict

NOSEQ = '*'


def _shift(value, delta):
    return value + delta if value is not None and value != -1 else value


class TagResult:
    """Tags detected for one molecule (one cDNA segment of a read)."""

    __slots__ = ('read_id', 'spatial_id', 'raw_spatial_id', 'molecule_id', 'bc_score', 'umi_good',
                 'strand', 'whitelist_index', 'polyT', 'primer', 'linker_start', 'linker_end', 'tso5')

    def __init__(self, read_id, spatial_id=NOSEQ, molecule_id=NOSEQ, bc_score=-1, umi_good=False,
                 strand='.', polyT=-1, primer=-1, linker_start=-1, linker_end=-1, tso5=-1,
                 raw_spatial_id=NOSEQ, whitelist_index=-1):
        self.read_id = read_id
        self.spatial_id = spatial_id or NOSEQ
        self.raw_spatial_id = raw_spatial_id or NOSEQ
        self.molecule_id = molecule_id or NOSEQ
        self.bc_score = bc_score
        self.umi_good = umi_good
        self.strand = strand
        self.whitelist_index = whitelist_index
        self.polyT = polyT
        self.primer = primer
        self.linker_start = linker_start
        self.linker_end = linker_end
        self.tso5 = tso5

    def is_valid(self):
        return self.spatial_id != NOSEQ

    def has_umi(self):
        return self.molecule_id != NOSEQ

    def set_strand(self, strand):
        self.strand = strand

    def update_coordinates(self, delta):
        self.polyT = _shift(self.polyT, delta)
        self.primer = _shift(self.primer, delta)
        self.linker_start = _shift(self.linker_start, delta)
        self.linker_end = _shift(self.linker_end, delta)
        self.tso5 = _shift(self.tso5, delta)

    def more_informative_than(self, other):
        if self.bc_score != other.bc_score:
            return self.bc_score > other.bc_score
        if self.linker_start != other.linker_start:
            return self.linker_start > other.linker_start
        if self.primer != other.primer:
            return self.primer > other.primer
        return self.polyT > other.polyT

    def segment_start(self):
        """Start of the molecule's segment in the oriented read (keeps some adapter upstream)."""
        return max(0, self.primer - 25, self.polyT - 75)

    def segment_end(self, seq_len):
        """End of the molecule's segment: just past the next TSO, or the end of the read."""
        if self.tso5 == -1:
            return seq_len
        return min(seq_len, self.tso5 + 25)

    def attributes(self):
        attr = []
        if self.polyT != -1:
            attr.append('PolyT detected')
        if self.primer != -1:
            attr.append('Primer detected')
        if self.linker_start != -1:
            attr.append('Linker detected')
        if self.tso5 != -1:
            attr.append('TSO detected')
        return attr

    def to_tsv(self):
        return '%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d' % (
            self.read_id, self.spatial_id, self.molecule_id, self.bc_score, self.umi_good, self.strand,
            self.polyT, self.primer, self.linker_start, self.linker_end, self.tso5)

    @staticmethod
    def header(spatial_label='spatial_id', molecule_label='molecule_id'):
        return ('read_id\t%s\t%s\tBC_score\tvalid_UMI\tstrand\t'
                'polyT_start\tprimer_end\tlinker_start\tlinker_end\tTSO5' % (spatial_label, molecule_label))


class SplitReadResult:
    """All molecules detected in one (possibly concatenated) read."""

    def __init__(self, read_id):
        self.read_id = read_id
        self.patterns = []

    def append(self, result):
        self.patterns.append(result)

    def empty(self):
        return not self.patterns

    def filter(self):
        """Keep only barcoded molecules; if none, keep the first result for reporting."""
        if not self.patterns:
            return
        valid = [r for r in self.patterns if r.is_valid()]
        self.patterns = valid if valid else [self.patterns[0]]


class TagStats:
    """Counters accumulated per worker and merged in the main process."""

    def __init__(self):
        self.counts = defaultdict(int)

    def add(self, key, value=1):
        self.counts[key] += value

    def add_result(self, result):
        self.counts['Molecules reported'] += 1
        for a in result.attributes():
            self.counts[a] += 1
        if result.is_valid():
            self.counts['Barcode detected'] += 1
        if result.has_umi():
            self.counts['UMI detected'] += 1

    def merge(self, other_counts):
        for k, v in other_counts.items():
            self.counts[k] += v

    def lines(self):
        return ['%s\t%d' % (k, self.counts[k]) for k in sorted(self.counts)]
