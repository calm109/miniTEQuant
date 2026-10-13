"""
Spatial-ID -> spot mapping (--barcode2spot).

Maps each spatial_id (e.g. Stereo-seq CID) to a spot: a segmented cell (cellbin),
a bin, a region, ... When a mapping is given, reads are tagged with the spot as
cell barcode so that `main.py quantify --st_mode` counts per spot instead of per CID.

File format (plain text or .gz), same syntax as spl-IsoQuant --barcode2spot:
    FILE                    spatial_id in column 0, spot in column 1
    FILE:ID_COL:SPOT_COL    0-based columns, e.g. cellbin.tsv:0:3
Columns are tab separated (comma or whitespace separated if the line has no tab).
Lines whose spatial_id is not a nucleotide sequence (e.g. a header) are skipped.
"""

import os
import re

from .seq_utils import open_text

_NUCLEOTIDES = set('ACGTN')
_COLUMNS_RE = re.compile(r'^(.+):(\d+):(\d+)$')


def parse_spec(spec):
    """Split 'FILE[:ID_COL:SPOT_COL]' into (path, id_col, spot_col)."""
    match = _COLUMNS_RE.match(spec)
    if match and not os.path.exists(spec):
        path, id_col, spot_col = match.group(1), int(match.group(2)), int(match.group(3))
    else:
        path, id_col, spot_col = spec, 0, 1
    if id_col == spot_col:
        raise ValueError('--barcode2spot: spatial ID and spot columns must differ: %s' % spec)
    if not os.path.exists(path):
        raise FileNotFoundError('--barcode2spot file not found: %s' % path)
    return path, id_col, spot_col


def _split(line):
    if '\t' in line:
        return [f.strip() for f in line.split('\t')]
    if ',' in line:
        return [f.strip() for f in line.split(',')]
    return line.split()


def sanitize_spot_name(name, name_sep):
    """
    Make a spot name safe as a miniTEQuant cell barcode: no whitespace, no name_sep
    (read names are split with rsplit(name_sep, 2)) and no trailing '-<digits>'
    (stripped as a 10x GEM suffix).
    """
    clean = re.sub(r'\s+', '.', name)
    if name_sep:
        clean = clean.replace(name_sep, '.')
    return re.sub(r'-(\d+)$', r'.\1', clean)


def load_spot_map(spec, name_sep, keep_ids=None):
    """
    Load the mapping as {spatial_id: spot}.

    keep_ids : optional set of spatial_ids to keep (the detected ones), which bounds
               memory for chip-wide mapping files.
    Returns (mapping, info) where info holds counters and renamed spot names.
    """
    path, id_col, spot_col = parse_spec(spec)
    need = max(id_col, spot_col) + 1
    mapping = {}
    spot_names = {}      # original -> sanitized, also interns spot strings
    info = {'lines': 0, 'skipped': 0, 'conflicts': 0, 'renamed': {}}
    with open_text(path) as handle:
        for line in handle:
            fields = _split(line.rstrip('\r\n'))
            if not fields or not fields[0]:
                continue
            info['lines'] += 1
            if len(fields) < need:
                info['skipped'] += 1
                continue
            spatial_id = fields[id_col].upper()
            spot = fields[spot_col]
            if not spatial_id or not set(spatial_id) <= _NUCLEOTIDES or not spot:
                info['skipped'] += 1      # header or malformed line
                continue
            if keep_ids is not None and spatial_id not in keep_ids:
                continue
            clean = spot_names.get(spot)
            if clean is None:
                clean = sanitize_spot_name(spot, name_sep)
                spot_names[spot] = clean
                if clean != spot:
                    info['renamed'][spot] = clean
            previous = mapping.get(spatial_id)
            if previous is None:
                mapping[spatial_id] = clean
            elif previous != clean:
                info['conflicts'] += 1    # keep the first assignment
    info['path'] = path
    info['columns'] = (id_col, spot_col)
    return mapping, info
