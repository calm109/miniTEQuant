"""
Spatial-barcode whitelist loading.

Accepted format (plain text or .gz), whitespace or comma separated:
    BARCODE                 one barcode per line
    BARCODE  X  Y           barcode with its spatial coordinates (e.g. Stereo-seq CID x y)
A header line is skipped automatically. All barcodes must have the same length (<= 32 nt).
"""

import numpy as np

from .seq_utils import encode_fixed_length, open_text

_NUCLEOTIDES = set('ACGTN')


class Whitelist:
    def __init__(self, barcodes, seq_len, x=None, y=None):
        self.barcodes = barcodes          # uint64 array, 2-bit encoded
        self.seq_len = seq_len
        self.x = x                        # int64 array or None
        self.y = y

    def __len__(self):
        return len(self.barcodes)

    @property
    def has_coordinates(self):
        return self.x is not None


def _split(line):
    return line.replace(',', ' ').split()


def load_whitelist(path, chunk_size=5000000):
    encoded_chunks, x_chunks, y_chunks = [], [], []
    barcodes, xs, ys = [], [], []
    seq_len = None
    has_coords = None
    skipped = 0

    def flush():
        encoded_chunks.append(encode_fixed_length(barcodes, seq_len))
        if has_coords:
            x_chunks.append(np.array(xs, dtype=np.int64))
            y_chunks.append(np.array(ys, dtype=np.int64))
        del barcodes[:], xs[:], ys[:]

    with open_text(path) as handle:
        for line in handle:
            fields = _split(line)
            if not fields:
                continue
            barcode = fields[0].upper()
            if not set(barcode) <= _NUCLEOTIDES:
                if seq_len is None:
                    continue    # header line
                skipped += 1
                continue
            if seq_len is None:
                seq_len = len(barcode)
                if seq_len > 32:
                    raise ValueError('Whitelist barcodes longer than 32 nt are not supported: %s' % barcode)
                has_coords = len(fields) >= 3
            if len(barcode) != seq_len:
                skipped += 1
                continue
            barcodes.append(barcode)
            if has_coords:
                xs.append(int(float(fields[1])))
                ys.append(int(float(fields[2])))
            if len(barcodes) >= chunk_size:
                flush()

    if seq_len is None:
        raise ValueError('No barcodes found in whitelist %s' % path)
    if barcodes:
        flush()

    encoded = np.concatenate(encoded_chunks)
    x = np.concatenate(x_chunks) if has_coords else None
    y = np.concatenate(y_chunks) if has_coords else None
    if skipped:
        print('[sp_tags] Whitelist: skipped %d malformed lines' % skipped, flush=True)
    return Whitelist(encoded, seq_len, x, y)
