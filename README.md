# ``treeCl`` - Phylogenetic Tree Clustering

``treeCl`` is a python package for clustering gene families by
phylogenetic similarity. It takes a collection of alignments, infers their phylogenetic trees,
and clusters them based on a matrix of between-tree distances. Finally, it calculates a single representative tree for each cluster.

You can read the paper [here](http://arxiv.org/abs/1510.02356)

## Installation

#### Preparing dependencies

There are some C and C++ dependencies that need to be available before installing
- [boost](http://www.boost.org/)
- [bio++](http://biopp.univ-montp2.fr/wiki/index.php/Main_Page)
- [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)

As well as some python dependencies
- numpy
- cython
- autowrap

The easiest way to get the C/C++ libraries is to use a package manager --- for instance homebrew for Mac (and Linux) --- otherwise there are installation instructions at the links above.

The python dependencies can all be installed using pip; numpy and cython are also available from conda.

#### External dependencies

To be able to build trees, treeCl needs to call on some external software. The choices are RAxML, PhyML, FastTree or PLL (using [pllpy](https://github.com/kgori/pllpy)). If any of these are installed, available in your path, and keep the standard names they were installed with, they should work.

#### Installing ``treeCl``
All remaining dependencies will be installed automatically using pip

    pip install treeCl



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
