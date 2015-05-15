#!/usr/bin/env python
import time
import treeCl
time.sleep(15)
c = treeCl.Collection(input_dir='/homes/kgori/scratch/simtest4', file_format='phylip')
c.calc_pll_trees()
