"""
Local alignment helpers used to locate adapters (primer / linker / TSO) and
to pick the best whitelist barcode for a read.

Scoring follows spl-IsoQuant: match +1, mismatch -1, gap -1 (linear).
Backend, in order of preference:
    1. ssw-py  (same library spl-IsoQuant uses)
    2. numba-compiled Smith-Waterman
    3. pure-Python Smith-Waterman
"""

import numpy as np

try:
    from ssw import AlignmentMgr
    BACKEND = 'ssw'
except ImportError:
    AlignmentMgr = None
    try:
        from numba import njit
        BACKEND = 'numba'
    except (ImportError, SystemError):
        BACKEND = 'python'


def _sw_python(query, ref):
    m = len(query)
    prev_h = [0] * (m + 1)
    prev_qs = [0] * (m + 1)
    prev_rs = [0] * (m + 1)
    best, bqs, bqe, brs, bre = 0, -1, -1, -1, -1
    for j, rc in enumerate(ref):
        cur_h = [0] * (m + 1)
        cur_qs = [0] * (m + 1)
        cur_rs = [0] * (m + 1)
        for i in range(1, m + 1):
            h, qs, rs = 0, 0, 0
            diag = prev_h[i - 1] + (1 if query[i - 1] == rc else -1)
            if diag > h:
                h = diag
                if prev_h[i - 1] > 0:
                    qs, rs = prev_qs[i - 1], prev_rs[i - 1]
                else:
                    qs, rs = i - 1, j
            up = prev_h[i] - 1
            if up > h:
                h, qs, rs = up, prev_qs[i], prev_rs[i]
            left = cur_h[i - 1] - 1
            if left > h:
                h, qs, rs = left, cur_qs[i - 1], cur_rs[i - 1]
            cur_h[i], cur_qs[i], cur_rs[i] = h, qs, rs
            if h > best:
                best, bqs, bqe, brs, bre = h, qs, i - 1, rs, j
        prev_h, prev_qs, prev_rs = cur_h, cur_qs, cur_rs
    return best, bqs, bqe, brs, bre


if BACKEND == 'numba':
    @njit(cache=True)
    def _sw_numba(query, ref):
        m = query.shape[0]
        prev_h = np.zeros(m + 1, np.int32)
        prev_qs = np.zeros(m + 1, np.int32)
        prev_rs = np.zeros(m + 1, np.int32)
        cur_h = np.zeros(m + 1, np.int32)
        cur_qs = np.zeros(m + 1, np.int32)
        cur_rs = np.zeros(m + 1, np.int32)
        best, bqs, bqe, brs, bre = 0, -1, -1, -1, -1
        for j in range(ref.shape[0]):
            cur_h[0] = 0
            for i in range(1, m + 1):
                h, qs, rs = 0, 0, 0
                diag = prev_h[i - 1] + (1 if query[i - 1] == ref[j] else -1)
                if diag > h:
                    h = diag
                    if prev_h[i - 1] > 0:
                        qs, rs = prev_qs[i - 1], prev_rs[i - 1]
                    else:
                        qs, rs = i - 1, j
                up = prev_h[i] - 1
                if up > h:
                    h, qs, rs = up, prev_qs[i], prev_rs[i]
                left = cur_h[i - 1] - 1
                if left > h:
                    h, qs, rs = left, cur_qs[i - 1], cur_rs[i - 1]
                cur_h[i], cur_qs[i], cur_rs[i] = h, qs, rs
                if h > best:
                    best, bqs, bqe, brs, bre = h, qs, i - 1, rs, j
            prev_h, cur_h = cur_h, prev_h
            prev_qs, cur_qs = cur_qs, prev_qs
            prev_rs, cur_rs = cur_rs, prev_rs
        return best, bqs, bqe, brs, bre


def local_align(query, ref):
    """
    Local alignment of query against ref.

    Returns (score, query_start, query_end, ref_start, ref_end); ends are inclusive.
    """
    if not query or not ref:
        return 0, -1, -1, -1, -1
    if BACKEND == 'ssw':
        mgr = AlignmentMgr(match_score=1, mismatch_penalty=1)
        mgr.set_read(query)
        mgr.set_reference(ref)
        a = mgr.align(gap_open=1, gap_extension=1)
        return a.optimal_score, a.read_start, a.read_end, a.reference_start, a.reference_end
    if BACKEND == 'numba':
        return _sw_numba(np.frombuffer(query.encode('ascii'), np.uint8),
                         np.frombuffer(ref.encode('ascii'), np.uint8))
    return _sw_python(query, ref)


# ---------------------------------------------------------------------------
# Pattern / barcode search (ported from spl-IsoQuant barcode_calling/common.py)
# ---------------------------------------------------------------------------

def align_pattern(sequence, start, end, pattern, min_score=0):
    """
    Align pattern to sequence[start:end].

    Returns (seq_start, seq_end, pattern_start, pattern_end, score) in read
    coordinates, or a tuple of Nones if the score is below min_score.
    """
    score, q_start, q_end, r_start, r_end = local_align(pattern, sequence[start:end])
    if score < min_score or q_start < 0:
        return None, None, None, None, None
    return start + r_start, start + r_end, q_start, q_end, score


def _best_window_alignment(sequence, start, end, kmer_size, pattern, positions, min_score, first_only):
    best = (None, None, None, None, 0)
    for match_position in positions:
        window_start = max(start, start + match_position - len(pattern) + kmer_size)
        window_end = min(end, start + match_position + len(pattern) + 1)
        alignment = align_pattern(sequence, window_start, window_end, pattern, min_score)
        if alignment[4] is None:
            continue
        if first_only:
            return alignment
        if alignment[4] > best[4]:
            best = alignment
    return best


def _finalize_positions(sequence, pattern, alignment, start_delta, end_delta, clip):
    start_pos, end_pos, pattern_start, pattern_end, _ = alignment
    if start_pos is None:
        return None, None
    if 0 <= start_delta < pattern_start:
        return None, None
    leftover_bases = len(pattern) - pattern_end - 1
    if 0 <= end_delta < leftover_bases:
        return None, None
    # Extend to the full pattern length, as if the unaligned pattern ends were present.
    if clip:
        return max(0, start_pos - pattern_start), min(end_pos + leftover_bases, len(sequence) - 1)
    return start_pos - pattern_start, end_pos + leftover_bases


def detect_exact_positions(sequence, start, end, kmer_size, pattern, positions,
                           min_score=0, start_delta=-1, end_delta=-1):
    """
    Locate pattern in sequence[start:end] using k-mer hit positions (relative to start).

    Takes the best-scoring alignment among all k-mer seeded windows.
    Returns (pattern_start, pattern_end) in read coordinates, or (None, None).
    """
    alignment = _best_window_alignment(sequence, start, end, kmer_size, pattern, positions,
                                       min_score, first_only=False)
    return _finalize_positions(sequence, pattern, alignment, start_delta, end_delta, clip=True)


def detect_first_exact_positions(sequence, start, end, kmer_size, pattern, positions,
                                 min_score=0, start_delta=-1, end_delta=-1):
    """Same as detect_exact_positions, but takes the first (leftmost) passing alignment."""
    alignment = _best_window_alignment(sequence, start, end, kmer_size, pattern, positions,
                                       min_score, first_only=True)
    return _finalize_positions(sequence, pattern, alignment, start_delta, end_delta, clip=False)


def find_best_barcode(candidates, read_barcode, min_score=10, score_diff=0, sufficient_score=0):
    """
    Pick the whitelist barcode that aligns best to the barcode region of a read.

    candidates: list of (barcode, shared_kmers, whitelist_index) from BarcodeIndex.query.
    Returns (barcode, score, start, end, whitelist_index) or (None, 0, 0, 0, -1).
    """
    best_barcode, best_index = None, -1
    best_score, best_start, best_end = 0, 0, 0
    second_best_score = 0

    for barcode, _, wl_index in candidates:
        score, q_start, q_end, r_start, r_end = local_align(barcode, read_barcode)
        if score < min_score:
            continue

        if score > best_score or (score == best_score and r_start - q_start < best_start):
            second_best_score = best_score
            best_barcode, best_index, best_score = barcode, wl_index, score
            best_start = r_start - q_start
            best_end = r_end + (len(barcode) - q_end)

        if score > sufficient_score > 0:
            # accept the first sufficiently good alignment
            break

    if best_score - second_best_score < score_diff:
        return None, 0, 0, 0, -1
    return best_barcode, best_score, best_start, best_end, best_index
