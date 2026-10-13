#!/usr/bin/env python
"""Command-line entry for post-alignment UMI deduplication (works when called as a script path)."""

import os
import sys

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from sp_tags.umi_dedup import main
else:
    from .umi_dedup import main

if __name__ == '__main__':
    main()
