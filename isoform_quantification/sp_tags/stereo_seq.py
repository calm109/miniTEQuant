"""
Stereo-seq tag detection for long reads (PacBio / ONT).

Read structure, in the orientation where the capture oligo comes first:

    PC1 primer | CID (25 nt) | linker | MID (10 nt) | polyT | cDNA (antisense) | TSO
      ...next molecule of a concatenated read: PC1 primer | CID | linker | MID | polyT | ...

Concatenated reads are split into molecules by scanning for successive
linker/polyT motifs; the TSO marks where one molecule ends.

Ported from spl-IsoQuant barcode_calling/callers/stereo.py
(StereoBarcodeDetector and StereoSplittingBarcodeDetector).
"""

from .align import detect_exact_positions, detect_first_exact_positions, find_best_barcode
from .base import SplitReadResult, TagResult
from .kmer_index import PatternIndex
from .seq_utils import find_polyt_start, reverse_complement


def _check_adapter(name, seq, kmer_size):
    seq = seq.upper()
    if not seq or set(seq) - set('ACGT'):
        raise ValueError('%s sequence must contain only A/C/G/T: %s' % (name, seq))
    if len(seq) < kmer_size + 3:
        raise ValueError('%s sequence is too short (%d nt, need >= %d): %s'
                         % (name, len(seq), kmer_size + 3, seq))
    return seq


class StereoSeqDetector:
    """Detects CID (spatial_id) and MID (molecule_id) in Stereo-seq long reads."""

    SPATIAL_LABEL = 'CID'
    MOLECULE_LABEL = 'MID'

    TSO5 = "CCCGCCTCTCAGTACGTCAGCAG"
    LINKER = "TTGTCTTCCTAAGAC"
    PC1_PRIMER = "CTTCCGATCTATGGCGACCTTATCAG"
    BC_LENGTH = 25
    UMI_LEN = 10
    TERMINAL_MATCH_DELTA = 3
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, barcode_index=None, split=True, min_score=None, primer=None, linker=None, tso=None):
        """
        barcode_index : BarcodeIndex over the CID whitelist, or None to report raw
                        (uncorrected) CIDs
        split         : True to split concatenated reads into molecules
        min_score     : minimal alignment score to accept a whitelist CID
                        (default: CID length - 3)
        primer, linker, tso : adapter sequences in capture-oligo-first orientation;
                        None uses PC1_PRIMER / LINKER / TSO5
        """
        self.primer = _check_adapter('primer', primer or self.PC1_PRIMER, 6)
        self.linker = _check_adapter('linker', linker or self.LINKER, 7)
        self.tso = _check_adapter('TSO', tso or self.TSO5, 8)

        # Thresholds tuned for the default adapters; capped at length - 3 so shorter
        # custom adapters stay detectable. Unchanged for the defaults.
        self.primer_min_score = min(12, len(self.primer) - 3)
        self.linker_min_score = min(12, len(self.linker) - 3)
        self.tso_min_score = min(15, len(self.tso) - 3)

        self.split = split
        self.barcode_index = barcode_index
        self.bc_length = barcode_index.seq_len if barcode_index is not None else self.BC_LENGTH
        self.min_score = min_score if min_score is not None else self.bc_length - 3

        # Splitting mode uses slightly stricter UMI / looser linker settings, as in spl-IsoQuant.
        self.umi_len_delta = 3 if split else 4
        self.relaxed_linker_min_score = min(10 if split else 12, len(self.linker) - 3)

        self.primer_index = PatternIndex(self.primer, 6)
        self.linker_index = PatternIndex(self.linker, 5)
        self.strict_linker_index = PatternIndex(self.linker, 7 if split else 6)
        self.tso5_index = PatternIndex(self.tso, 8)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def detect(self, read_id, sequence):
        """Return the list of TagResult for a read (several if the read is concatenated)."""
        if self.split:
            return self._detect_split(read_id, sequence).patterns
        return [self._detect_single(read_id, sequence)]

    # ------------------------------------------------------------------
    # Read-level search
    # ------------------------------------------------------------------

    def _detect_single(self, read_id, sequence):
        result = self._find_fwd(read_id, sequence)
        if result.polyT != -1:
            result.set_strand('+')
        if result.is_valid():
            return result

        rev_result = self._find_fwd(read_id, reverse_complement(sequence))
        if rev_result.polyT != -1:
            rev_result.set_strand('-')
        if rev_result.is_valid():
            return rev_result

        return result if result.more_informative_than(rev_result) else rev_result

    def _detect_split(self, read_id, sequence):
        read_result = SplitReadResult(read_id)
        r = None
        for strand, oriented in (('+', sequence), ('-', reverse_complement(sequence))):
            r = self._find_fwd(read_id, oriented)
            prev_start = 0
            while r.polyT != -1:
                r.set_strand(strand)
                read_result.append(r)
                if r.tso5 != -1:
                    current_start = r.tso5 + 15
                else:
                    current_start = r.polyT + 100
                # always make a step
                current_start = max(prev_start + 150, current_start)
                prev_start = current_start
                if len(oriented) - current_start < 50:
                    break
                r = self._find_fwd(read_id, oriented[current_start:])
                r.update_coordinates(current_start)

        if read_result.empty():
            read_result.append(r)
        read_result.filter()
        return read_result

    # ------------------------------------------------------------------
    # Single-molecule search in forward orientation
    # ------------------------------------------------------------------

    def _find_tso(self, sequence, polyt_start):
        search_start = polyt_start + 1
        positions = self.tso5_index.find(sequence, search_start, len(sequence))
        tso_start, _ = detect_first_exact_positions(sequence, search_start, len(sequence),
                                                    self.tso5_index.k, self.tso, positions, min_score=self.tso_min_score,
                                                    start_delta=self.TERMINAL_MATCH_DELTA,
                                                    end_delta=self.TERMINAL_MATCH_DELTA)
        return tso_start

    def _find_fwd(self, read_id, sequence):
        polyt_start = find_polyt_start(sequence)

        linker_start, linker_end = None, None
        tso5_start = None
        if polyt_start != -1:
            # use relaxed parameters if polyT is found
            positions = self.linker_index.find(sequence, 0, polyt_start + 1)
            linker_start, linker_end = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                              self.linker_index.k, self.linker, positions,
                                                              min_score=self.relaxed_linker_min_score,
                                                              start_delta=self.TERMINAL_MATCH_DELTA,
                                                              end_delta=self.TERMINAL_MATCH_DELTA)
            if self.split:
                tso5_start = self._find_tso(sequence, polyt_start)

        if linker_start is None:
            # polyT not found, or no linker upstream of it: look for the linker in the entire read
            detect = detect_first_exact_positions if self.split else detect_exact_positions
            positions = self.strict_linker_index.find(sequence)
            linker_start, linker_end = detect(sequence, 0, len(sequence),
                                              self.strict_linker_index.k, self.linker, positions, min_score=self.linker_min_score,
                                              start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                              end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if linker_start is None:
            return TagResult(read_id, polyT=polyt_start)

        if polyt_start == -1 or (self.split and polyt_start < linker_start):
            # once the linker is found, look for a short polyT right after the MID
            presumable_polyt_start = linker_end + self.UMI_LEN
            search_start = max(0, presumable_polyt_start - 4)
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start
            if self.split:
                tso5_start = self._find_tso(sequence, polyt_start)

        if tso5_start is not None:
            # another linker between polyT and TSO means the TSO belongs to a later molecule
            positions = self.strict_linker_index.find(sequence, polyt_start + 1, tso5_start)
            new_linker_start, _ = detect_exact_positions(sequence, polyt_start + 1, tso5_start,
                                                         self.strict_linker_index.k, self.linker, positions,
                                                         min_score=self.linker_min_score,
                                                         start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                         end_delta=self.STRICT_TERMINAL_MATCH_DELTA)
            if new_linker_start is not None and new_linker_start - polyt_start > 100:
                tso5_start = new_linker_start - self.bc_length - len(self.primer) - len(self.tso)
        else:
            tso5_start = -1

        positions = self.primer_index.find(sequence, 0, linker_start)
        primer_start, primer_end = detect_exact_positions(sequence, 0, linker_start,
                                                          self.primer_index.k, self.primer, positions,
                                                          min_score=self.primer_min_score, end_delta=self.TERMINAL_MATCH_DELTA)
        if primer_start is None:
            primer_end = linker_start - self.bc_length - 1

        layout = dict(polyT=polyt_start, linker_start=linker_start, linker_end=linker_end, tso5=tso5_start)
        if primer_end < 0:
            return TagResult(read_id, primer=-1, **layout)
        layout['primer'] = primer_end

        barcode_start = primer_end + 1
        barcode_end = linker_start - 1
        if abs(barcode_end - barcode_start - self.bc_length) > 10:
            return TagResult(read_id, **layout)

        raw_barcode = sequence[barcode_start:barcode_end + 1]
        if self.barcode_index is None:
            barcode, bc_score, wl_index = raw_barcode, 0, -1
        else:
            candidates = self.barcode_index.query(raw_barcode, max_hits=10, min_kmers=2)
            barcode, bc_score, _, _, wl_index = find_best_barcode(candidates, raw_barcode,
                                                                  min_score=self.min_score,
                                                                  sufficient_score=self.bc_length - 1)
            if barcode is None:
                return TagResult(read_id, raw_spatial_id=raw_barcode, **layout)

        umi_start = linker_end + 1
        umi_end = polyt_start - 1
        if umi_end - umi_start <= 5:
            umi_end = umi_start + self.UMI_LEN - 1
        umi = None
        good_umi = False
        if umi_start + 2 * self.UMI_LEN > umi_end > umi_start:
            umi = sequence[umi_start:umi_end + 1]
            if len(umi) < self.UMI_LEN and polyt_start != -1:
                # polyT start snaps to the first 'TTTT', swallowing MID bases that end in T
                full_umi = sequence[umi_start:umi_start + self.UMI_LEN]
                if not full_umi[len(umi):].strip('T'):
                    umi = full_umi
            good_umi = abs(len(umi) - self.UMI_LEN) <= self.umi_len_delta

        return TagResult(read_id, barcode, umi, bc_score, good_umi, raw_spatial_id=raw_barcode,
                         whitelist_index=wl_index, **layout)
