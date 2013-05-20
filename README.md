# ``TreeCl`` - Phylogenetic Tree Clustering

TreeCl is a python package for clustering gene families by
phylogenetic similarity. It takes a collection of alignments, infers their phylogenetic trees, 
and clusters them based on a matrix of between-tree distances. Finally, it calculates a single representative tree for each cluster.

The purpose of this is to establish whether there is any underlying structure
to the data.

## Installation 

Clone the repo with submodules using
 ```git clone --recursive git@github.com:kgori/treeCl.git``` 
 and add it to your $PYTHONPATH

## Dependencies

#### Python:
- [numpy](http://www.numpy.org "NumPy") (v1.6.2) 
- [scipy](http://www.scipy.org "SciPy") (v0.11.0)
- [dendropy](http://pythonhosted.org/DendroPy/ "DendroPy is a Python library for phylogenetic computing.") (v3.12.0)
- [scikit-learn](http://scikit-learn.org/stable/ "Machine learning in Python") (v0.12.1)

The easiest way to install the dependencies is using pip. If you don't have pip,
it can be installed by typing ```easy_install pip``` in a shell.
Then the above packages can be installed by running this command:
    
    pip install numpy scipy dendropy scikit-learn

#### External:
- [ALF](http://darwin-services.inf.ethz.ch/DarwinServices/ALF.html#service0 "ALF: simulating genome evolution") 
- [Darwin](http://www.cbrg.ethz.ch/darwin "Data Analysis and Retrieval With Indexed Nucleotide/peptide sequences")
- [PhyML](https://code.google.com/p/phyml/ "Phylogenetic estimation using Maximum Likelihood")
- TreeCollection

#### Other:
- [GTP](http://dl.acm.org/citation.cfm?id=1916603 "ACM digital library") - a java program for calculating geodesic distances - see [A Fast Algorithm for Computing Geodesic Distances in Tree Space](https://cs.uwaterloo.ca/~m2owen/pub/poly_geodesic.pdf "Owen and Provan, 2011") 

## Example Analysis
```
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Clustering, Partition

c = Collection(input_dir='input_dir', file_format='phylip', datatype='protein') # add compression='gz' or 'bz2' if sequence alignments are compressed (zip not supported yet)
c.calc_NJ_trees() #add verbosity=1 or higher to get progress messages
dm = c.distance_matrix('euc')
cl = Clustering(dm)
p = cl.hierarchical(4, 'single') # should give fairly inaccurate clustering
true = Partition(tuple([1]*15+[2]*15+[3]*15+[4]*15))
sc = Scorer(c.records)
score = sc.score(p)
print score
```



