# standard library
import os
import re

# third party
import numpy as np

TMPDIR = os.getenv('TMPDIR', '/tmp')
EPS = 1e-8 # a small number
SORT_KEY = lambda item: tuple((int(num) if num else alpha) for (num, alpha) in
                              re.findall(r'(\d+)|(\D+)', item))
ANALYSES = ['tlr', 'lr', 'l', 'r', 'ml', 'full', 'nj', 'bionj', 'bionj+', 'lk']
DNA_ACGT_THRESHOLD = 0.75 # proportion of ACGT in sequence to call it as DNA
POSINF = float('inf')     # positive infinity
NEGINF = -POSINF          # negative infinity