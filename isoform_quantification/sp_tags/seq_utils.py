"""
Sequence utilities for sp_tags: polyT detection, reverse complement,
2-bit encoding of barcodes and streaming readers for FASTA / FASTQ / BAM.

polyT detection and 2-bit encoding follow spl-IsoQuant
(isoquant_lib/barcode_calling/common.py).
"""

import gzip
import os

import numpy as np


# ---------------------------------------------------------------------------
# Basic sequence operations
# ---------------------------------------------------------------------------

_COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def reverse_complement(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def find_polyt_start(seq, window_size=16, polya_fraction=0.75):
    """
    Return the start of the first polyT stretch in seq, or -1.

    A sliding window of window_size is moved along seq until it contains at
    least window_size * polya_fraction T's; the start is then snapped to the
    first 'TTTT' inside the window.
    """
    polya_count = int(window_size * polya_fraction)

    if len(seq) < window_size:
        return -1
    i = 0
    t_count = seq[0:window_size].count('T')
    while i < len(seq) - window_size:
        if t_count >= polya_count:
            break
        first_base_t = seq[i] == 'T'
        new_base_t = i + window_size < len(seq) and seq[i + window_size] == 'T'
        if first_base_t and not new_base_t:
            t_count -= 1
        elif not first_base_t and new_base_t:
            t_count += 1
        i += 1

    if i >= len(seq) - window_size:
        return -1

    return i + max(0, seq[i:].find('TTTT'))


# ---------------------------------------------------------------------------
# 2-bit encoding (A=0, C=1, T=2, G=3 via (ord(c) & 6) >> 1; N is read as G)
# ---------------------------------------------------------------------------

BIN2NUCL = "ACTG"


def str_to_2bit(seq):
    value = 0
    for c in seq:
        value = (value << 2) | ((ord(c) & 6) >> 1)
    return value


def bit_to_str(value, seq_len):
    value = int(value)
    chars = []
    for i in range(seq_len):
        chars.append(BIN2NUCL[(value >> ((seq_len - i - 1) * 2)) & 3])
    return ''.join(chars)


def encode_fixed_length(barcodes, seq_len):
    """Encode a list of equal-length strings into a uint64 array (seq_len <= 32)."""
    n = len(barcodes)
    if n == 0:
        return np.zeros(0, dtype=np.uint64)
    raw = np.frombuffer(''.join(barcodes).encode('ascii'), dtype=np.uint8).reshape(n, seq_len)
    codes = ((raw & 6) >> 1).astype(np.uint64)
    result = np.zeros(n, dtype=np.uint64)
    for pos in range(seq_len):
        result |= codes[:, pos] << np.uint64((seq_len - pos - 1) * 2)
    return result


def kmer_codes(seq, k):
    """2-bit codes of all k-mers of seq (sliding window), as a uint64 array."""
    n = len(seq) - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64)
    codes = ((np.frombuffer(seq.encode('ascii'), dtype=np.uint8) & 6) >> 1).astype(np.uint64)
    result = np.zeros(n, dtype=np.uint64)
    for i in range(k):
        result |= codes[i:i + n] << np.uint64((k - 1 - i) * 2)
    return result


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def _strip_compression(path):
    base = os.path.basename(path)
    lower = base.lower()
    for ext in ('.gz', '.gzip'):
        if lower.endswith(ext):
            return base[:-len(ext)], True
    return base, False


def detect_format(path):
    """Return 'fastq', 'fasta' or 'bam' from the file extension."""
    base, _ = _strip_compression(path)
    ext = os.path.splitext(base)[1].lower()
    if ext in ('.fq', '.fastq'):
        return 'fastq'
    if ext in ('.fa', '.fasta', '.fna'):
        return 'fasta'
    if ext in ('.bam', '.sam'):
        return 'bam'
    raise ValueError('Unknown read file format: %s (expected FASTQ, FASTA or BAM)' % path)


def open_text(path):
    if _strip_compression(path)[1]:
        return gzip.open(path, 'rt')
    return open(path, 'r')


def _iter_fastq(path):
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            if not header.strip():
                continue
            seq = handle.readline().rstrip()
            handle.readline()
            qual = handle.readline().rstrip()
            read_id = header[1:].split(None, 1)[0]
            yield read_id, seq.upper(), qual


def _iter_fasta(path):
    with open_text(path) as handle:
        read_id = None
        chunks = []
        for line in handle:
            if line.startswith('>'):
                if read_id is not None:
                    yield read_id, ''.join(chunks).upper(), None
                read_id = line[1:].split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if read_id is not None:
            yield read_id, ''.join(chunks).upper(), None


def _iter_bam(path):
    import pysam
    mode = 'r' if path.lower().endswith('.sam') else 'rb'
    with pysam.AlignmentFile(path, mode, check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            seq = r.query_sequence
            if not seq:
                continue
            quals = r.query_qualities
            qual = None if quals is None else ''.join(chr(q + 33) for q in quals)
            yield r.query_name, seq.upper(), qual


def iter_reads(path):
    """Yield (read_id, sequence, quality_or_None) from FASTQ / FASTA / BAM (optionally gzipped)."""
    fmt = detect_format(path)
    if fmt == 'fastq':
        return _iter_fastq(path)
    if fmt == 'fasta':
        return _iter_fasta(path)
    return _iter_bam(path)
