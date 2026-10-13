"""
sp_tags: spatial-barcode and UMI processing for high-resolution spatial long-read data.

Platform identifiers are mapped onto two unified fields:
    spatial_id  : Stereo-seq CID, Visium HD barcode, ...
    molecule_id : Stereo-seq MID, 10x UMI, ...

Run:  python isoform_quantification/sp_tags/run_sp_tags.py -h
"""

from .base import TagResult, TagStats
from .kmer_index import BarcodeIndex
from .stereo_seq import StereoSeqDetector
from .umi import cluster_umis
from .whitelist import load_whitelist
