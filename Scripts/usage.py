import sys
sys.path.append('/home/malcolm/Documents/EBi/Git/')
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Clustering, Partition
from treeCl.em import emtrees

# Collected Usage:

c = Collection(input_dir='/home/malcolm/Documents/EBi/Data/malcolm_test_data', compression='gz', file_format='phylip', datatype='protein')
c.calc_NJ_trees() #add verbosity=1 or higher to get progress messages
dm = c.distance_matrix('euc')
# cl = Clustering(dm)
# p = cl.hierarchical(4, 'single') # should give fairly inaccurate clustering
# true = Partition(tuple([1]*15+[2]*15+[3]*15+[4]*15))
sc = Scorer(c.records,c.analysis)

# sc.score(true)
# sc.score(p)
e = emtrees(sc,4)

test = open('test1', 'w')

for i in range(1):
    e.assign_clusters()
    e.maximise()
    test.write(str(e.partition.partition_vector) + '\n')

