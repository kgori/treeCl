#!/usr/bin/env python
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Clustering, Partition
from treeCl.algorithms import emtrees
from multiprocessing import Pool

c = Collection(input_dir='/home/malcolm/Documents/EBi/Data/malcolm_test_data/',
               compression='gz', file_format='phylip', datatype='protein')
c.calc_NJ_trees()

scores = open('scores', 'w')
clusters = open('clusters', 'w')

def trees(method):
    e = emtrees(c, 4)
    e.maximise(method)
    scores.write(str(e.L) + '\n')
    clusters.write(e.partition.partition_vector + '\n')

pool = Pool(2)

scores.write('Liklihood:')
clusters.write('Liklihood:')

scores1 = map(trees, ['ml'] * 10)

scores.writei('Distance:')
clusters.write('Distance:')

scores2 = pool.map(trees, ['ml'] * 10)

print(scores1)
print(scores2)
