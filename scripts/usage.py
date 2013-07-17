#!/usr/bin/env python
import sys
sys.path.append('/home/malcolm/Documents/EBi/Git/')
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Clustering, Partition
from treeCl.algorithms import EMTrees
from analysis import Result

# Collection Usage:

c = Collection(input_dir='/home/malcolm/Documents/EBi/Data/test_data',
               compression='gz', file_format='phylip', datatype='protein')
c.calc_NJ_trees()
dm = c.distance_matrix('geo')
# cl = Clustering(dm)
# sd = cl.spectral_decomp('estimate', 'estimate')
# p = cl.spectral_cluster(4, sd)
# true = Partition(tuple([1]*15+[2]*15+[3]*15+[4]*15))
# sc = Scorer(c.records, c.analysis)

# # sc.score(true)
# score = sc.score(p)
# result = Result(score, p.get_membership(), sc.history())

# emtrees usage:

e = EMTrees(c, 4)
e.random_partition()
# e.maximise_test('dist')
# print e.L
