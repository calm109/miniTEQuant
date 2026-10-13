#!/usr/bin/env python
"""Command-line entry for sp_tags (works when called as a script path)."""

import os
import sys

if __package__ in (None, ''):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from sp_tags.pipeline import main
else:
    from .pipeline import main

if __name__ == '__main__':
    main()
