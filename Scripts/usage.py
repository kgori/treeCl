import sys
sys.path.append('/home/malcolm/Documents/EBi/Git/')
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Clustering, Partition
from treeCl.algorithms import emtrees

# Collection Usage:

c = Collection(input_dir='/home/malcolm/Documents/EBi/Data/easy_case', compression='gz', file_format='phylip', datatype='protein')
c.calc_NJ_trees() #add verbosity=1 or higher to get progress messages
dm = c.distance_matrix('euc')
cl = Clustering(dm)
p = cl.hierarchical(4, 'single') # should give fairly inaccurate clustering
true = Partition(tuple([1]*15+[2]*15+[3]*15+[4]*15))
sc = Scorer(c.records,c.analysis)

sc.score(true)
sc.score(p)

e = emtrees(c,4)
# e.maximise('ml')
# print e.L