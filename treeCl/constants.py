# standard library
import os
import re

# third party
import numpy as np

TMPDIR = os.getenv('TMPDIR', '/tmp')
EPS = 1e-8
MINUS_INF = -np.inf
SORT_KEY = lambda item: tuple((int(num) if num else alpha) for (num, alpha) in
                              re.findall(r'(\d+)|(\D+)', item))
ANALYSES = ['tlr', 'lr', 'l', 'r', 'ml', 'full', 'nj', 'bionj', 'bionj+', 'lk']
DNA_ACGT_THRESHOLD = 0.75
