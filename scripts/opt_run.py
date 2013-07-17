#!/usr/bin/env python

from collection import Collection
from optimiser import Optimiser
import sys
import tempfile
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

new_tmpdir = tempfile.mkdtemp(prefix='tmpwrap_mgp_', dir='/tmp/')

c = Collection(input_dir='/nfs/gns/homes/mgperry/treeCl_data/test_data/',
               compression='gz', file_format='phylip', datatype='protein',
               tmpdir=new_tmpdir)

o = Optimiser(4, c)
o.optimise(max_iter=500, nreassign = 10)

scores = open('output_' + id_generator(6), 'w+')

cl_dict = o.get_clusters(o.global_best_assignment)

scores.write('Likelihood score:' + str(o.global_best_score) + '\n')
scores.write('Clusters:')
scores.write('\n'.join([str(cl_dict[cl]) for cl in sorted(cl_dict.keys())]))
scores.write('\n\n')

o.Scorer.print_history(scores)
