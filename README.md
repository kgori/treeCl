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
- [biopython](http://www.biopython.org/‎ "Biopython") (v1.60) *optional - for k-medoids clustering only*

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
``` python
"""
Import some classes from treeCl
"""
from treeCl import Collection, Scorer, Clustering, Partition, DistanceMatrix


"""
Load your data. This should be a directory full of sequence alignments in fasta '*.fas'
or phylip '*.phy' formats.
"""
c = Collection(input_dir='input_dir', file_format='phylip', 
    datatype='protein') # add compression='gz' or 'bz2' 
                        # if sequence alignments are compressed 
                        # (zip not supported yet)

"""
Calculate some trees. Trees can be Neighbour Joining ('nj') or maximum 
likelihood ('ml'), or an intermediate algorithm where the topology is 
done by NJ, and the branch lengths are done by ML.
The intermediate algorithm is run with 
    c.calc_NJ_trees(analysis='lr') # 'lr' = 'lengths and rates'
"""
c.calc_NJ_trees() # add verbosity=1 or higher to get progress messages
# c.calc_ML_trees() # use maximum likelihood - slower, more accurate

"""
Now get the distance matrix of between-tree distances and set up a clustering object
"""
dm = c.distance_matrix('euc')
cl = Clustering(dm)


"""
Do some clustering. Syntax is
    cl.hierarchical(<number_of_clusters:int>, <method:str>)
        'method' is one of 'single', 'complete', 'average', 'ward'
    cl.kmedoids(<number_of_clusters:int>)
    cl.spectral_cluster(<number_of_clusters:int>, <decomp:Decomp>)
    cl.MDS_cluster(<number_of_clusters:int>, <decomp:Decomp>)

Spectral and MDS methods need a 'Decomp' object, which is an eigen decomposition
into eigenvalues and eigenvectors. These are obtained with:
    decomp = cl.spectral_decomp(-1, 'median') # The -1 and 'median' parameters seem to work best
    decomp = cl.MDS_decomp()
"""

p = cl.hierarchical(2, 'single') # example single-linkage hierarchical 
                                 # clustering into 2 groups

"""
Score the result via likelihood
"""
sc = Scorer(c.records, analysis='nj') # or 'ml' or 'lr'
score = sc.score(p)
print score

"""
Finally, get the sequences and trees for the two groups
"""

groups = p.get_membership()
sequences1 = sc.concatenate(groups[0])  # sequences1.write_fasta('filename') 
                                        # or .write_phylip('filename', 
                                        #         interleaved=True)
sequences2 = sc.concatenate(groups[1])

tree1 = sc.concats[groups[0]]
tree2 = sc.concats[groups[1]]

```



