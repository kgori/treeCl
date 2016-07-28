# -*- coding: utf-8 -*-
# standard library
import os
import re
import sys

TMPDIR = os.getenv('TMPDIR', '/tmp')
EPS = 1e-8  # a small number
SORT_KEY = lambda item: tuple((int(num) if num else alpha) for (num, alpha) in
                              re.findall(r'(\d+)|(\D+)', item))
ANALYSES = ['tlr', 'lr', 'l', 'r', 'ml', 'full', 'nj', 'bionj', 'bionj+', 'lk']
DNA_ACGT_THRESHOLD = 0.75  # proportion of ACGT in sequence to call it as DNA
POSINF = float('inf')  # positive infinity
NEGINF = -POSINF  # negative infinity
PARALLEL_PROFILE = os.getenv('IPYTHON_PROFILE', None)
PHYML_MEMORY_MULTIPLIER = 1.2  # Add %age wiggle room when calculating Phyml
# memory requirements on LSF
PHYML_MEMORY_SPARE = 256  # Allocate extra Mb of memory on LSF jobs
PHYML_MEMORY_MIN = 512
PHYML_MEMORY_STANDARD = 4096
RANDOM_SEED = int('CA55E77E', 16)
logo = """
═══════════ ╔═╗┬
┌┬┐┬─┐┌─┐┌─┐║  │
 │ ├┬┘├┤ ├┤ ╚═╝┴─┘
 ┴ ┴└─└─┘└─┘╭─────
┈┈┈┈┈┈┄┄┄┄┄─┤  ╭──
            ╰──┤
══════════════ ╰──
"""
ISPY3 = sys.version_info.major>=3