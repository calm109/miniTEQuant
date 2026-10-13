"""
UMI (molecule_id) error correction.

Long reads carry sequencing errors in the UMI, so copies of one molecule can
show slightly different UMIs. Within each spatial_id, UMIs are clustered
greedily: UMIs are visited from most to least frequent, and each one is merged
into the closest already accepted UMI within max_edit_distance, otherwise it is
accepted as a new molecule.

The greedy scheme follows spl-IsoQuant UMIFilter._construct_umi_dict. It is
applied per spatial_id (before alignment, so without gene information); for
Stereo-seq this is safe because one CID (~0.5 um spot) holds few molecules.
"""

try:
    import editdistance
    _edit_distance = editdistance.eval
except ImportError:
    _edit_distance = None


def bounded_edit_distance(a, b, max_dist):
    """Levenshtein distance of a and b, or max_dist + 1 if it exceeds max_dist."""
    if abs(len(a) - len(b)) > max_dist:
        return max_dist + 1
    if _edit_distance is not None:
        return _edit_distance(a, b)
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        row_min = i
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
            row_min = min(row_min, cur[j])
        if row_min > max_dist:
            return max_dist + 1
        prev = cur
    return prev[-1]


def cluster_umis(umi_counts, max_edit_distance=1, umi_length=0):
    """
    Cluster the UMIs observed for one spatial_id.

    umi_counts : dict {umi: number of reads}
    umi_length : expected UMI length; on equal counts, UMIs closer to it win (0 = ignore)

    Returns dict {umi: representative_umi}.
    """
    if max_edit_distance <= 0 or len(umi_counts) < 2:
        return {u: u for u in umi_counts}

    def sorting_key(umi):
        length_penalty = -abs(umi_length - len(umi)) if umi_length else 0
        return umi_counts[umi], length_penalty, umi

    trusted = []
    mapping = {}
    for umi in sorted(umi_counts, key=sorting_key, reverse=True):
        best_umi, best_dist = None, max_edit_distance + 1
        for candidate in trusted:
            dist = bounded_edit_distance(umi, candidate, max_edit_distance)
            if dist < best_dist:
                best_umi, best_dist = candidate, dist
                if dist <= 1:
                    break
        if best_umi is None:
            trusted.append(umi)
            mapping[umi] = umi
        else:
            mapping[umi] = best_umi
    return mapping
