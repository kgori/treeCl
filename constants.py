import numpy as np
import os
import re

TMPDIR = os.getenv('TMPDIR', '/tmp')
EPS = 1e-8
MINUS_INF = -np.inf
SORT_KEY = lambda item: tuple((int(num) if num else alpha) for (num, alpha) in
                              re.findall(r'(\d+)|(\D+)', item))
ANALYSES = ['tlr', 'lr', 'l', 'r', 'ml', 'full', 'nj', 'bionj', 'bionj+', 'lk']
