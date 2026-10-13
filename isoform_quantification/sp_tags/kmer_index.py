"""
K-mer indexes.

PatternIndex : finds k-mer seeds of a single adapter (primer / linker / TSO) in a read.
BarcodeIndex : k-mer index over a spatial-barcode whitelist, stored as numpy arrays
               so that it scales to Stereo-seq chips with tens of millions of CIDs and
               is shared copy-on-write by forked worker processes.
"""

import re

import numpy as np

from .seq_utils import bit_to_str, kmer_codes


class PatternIndex:
    """K-mer seeds of one adapter sequence, matched with a compiled regex."""

    def __init__(self, pattern, kmer_size):
        self.pattern = pattern
        self.k = kmer_size
        kmers = sorted({pattern[i:i + kmer_size] for i in range(len(pattern) - kmer_size + 1)})
        self._regex = re.compile('(?=(?:%s))' % '|'.join(kmers))

    def find(self, sequence, start=0, end=None):
        """Positions (relative to start) of k-mers of the pattern inside sequence[start:end]."""
        if end is None:
            end = len(sequence)
        start = max(0, start)
        if end - start < self.k:
            return []
        return [m.start() - start for m in self._regex.finditer(sequence, start, end)]


class BarcodeIndex:
    """
    K-mer index over 2-bit encoded whitelist barcodes.

    Two layouts:
      * dense  : a 4^k offset table (direct addressing), for very large whitelists
      * sparse : sorted unique k-mers + offsets (binary search), for smaller ones
    """

    DENSE_FILL_FACTOR = 0.25
    BUILD_CHUNK = 5000000

    def __init__(self, encoded_barcodes, seq_len, kmer_size=14):
        self.barcodes = np.asarray(encoded_barcodes, dtype=np.uint64)
        self.seq_len = seq_len
        self.k = kmer_size
        self.mask = np.uint64((1 << (2 * kmer_size)) - 1)
        self.kmers_per_barcode = seq_len - kmer_size + 1
        self.shifts = np.array([(seq_len - kmer_size - i) * 2 for i in range(self.kmers_per_barcode)],
                               dtype=np.uint64)
        n = len(self.barcodes)
        self.id_dtype = np.uint32 if n < (1 << 32) else np.uint64

        total_entries = n * self.kmers_per_barcode
        self.dense = total_entries > (4 ** kmer_size) * self.DENSE_FILL_FACTOR
        if self.dense:
            self._build_dense()
        else:
            self._build_sparse()

    def __len__(self):
        return len(self.barcodes)

    def _chunk_kmers(self, start, end):
        chunk = self.barcodes[start:end]
        return ((chunk[:, None] >> self.shifts[None, :]) & self.mask).ravel()

    def _build_sparse(self):
        n = len(self.barcodes)
        kmers = self._chunk_kmers(0, n)
        order = np.argsort(kmers, kind='stable')
        self.index = (order // self.kmers_per_barcode).astype(self.id_dtype)
        self.keys, starts = np.unique(kmers[order], return_index=True)
        self.offsets = np.append(starts, len(kmers)).astype(np.int64)

    def _build_dense(self):
        n = len(self.barcodes)
        total_kmers = 4 ** self.k
        counts = np.zeros(total_kmers, dtype=np.int64)
        for start in range(0, n, self.BUILD_CHUNK):
            kmers = self._chunk_kmers(start, min(n, start + self.BUILD_CHUNK)).astype(np.int64)
            counts += np.bincount(kmers, minlength=total_kmers)

        self.offsets = np.zeros(total_kmers + 1, dtype=np.int64)
        np.cumsum(counts, out=self.offsets[1:])
        del counts
        self.index = np.empty(int(self.offsets[-1]), dtype=self.id_dtype)

        write_ptr = self.offsets[:-1].copy()
        for start in range(0, n, self.BUILD_CHUNK):
            kmers = self._chunk_kmers(start, min(n, start + self.BUILD_CHUNK)).astype(np.int64)
            order = np.argsort(kmers, kind='stable')
            sorted_kmers = kmers[order]
            uniq, first, cnt = np.unique(sorted_kmers, return_index=True, return_counts=True)
            rank = np.arange(len(sorted_kmers), dtype=np.int64) - np.repeat(first, cnt)
            self.index[write_ptr[sorted_kmers] + rank] = start + order // self.kmers_per_barcode
            write_ptr[uniq] += cnt

    def _kmer_ranges(self, query_kmers):
        if self.dense:
            q = query_kmers.astype(np.int64)
            return self.offsets[q], self.offsets[q + 1]
        pos = np.searchsorted(self.keys, query_kmers)
        found = pos < len(self.keys)
        found[found] = self.keys[pos[found]] == query_kmers[found]
        pos = pos[found]
        return self.offsets[pos], self.offsets[pos + 1]

    def query(self, sequence, max_hits=0, min_kmers=1, hits_delta=1):
        """
        Find whitelist barcodes sharing k-mers with sequence.

        Returns a list of (barcode, shared_kmers, whitelist_index), best first,
        keeping only hits within hits_delta of the top shared k-mer count.
        """
        query_kmers = kmer_codes(sequence, self.k)
        if len(query_kmers) == 0 or len(self.barcodes) == 0:
            return []
        starts, ends = self._kmer_ranges(query_kmers)
        hits = [self.index[s:e] for s, e in zip(starts, ends) if e > s]
        if not hits:
            return []
        ids, counts = np.unique(np.concatenate(hits), return_counts=True)
        keep = counts >= max(min_kmers, counts.max() - hits_delta)
        ids, counts = ids[keep], counts[keep]
        order = np.argsort(-counts, kind='stable')
        if max_hits:
            order = order[:max_hits]
        return [(bit_to_str(self.barcodes[ids[i]], self.seq_len), int(counts[i]), int(ids[i]))
                for i in order]
