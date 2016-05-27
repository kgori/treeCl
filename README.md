# ``treeCl`` - Phylogenetic Tree Clustering

![Travis build status](https://travis-ci.org/kgori/treeCl.svg?branch=travis)

``treeCl`` is a python package for clustering gene families by
phylogenetic similarity. It takes a collection of alignments, infers their phylogenetic trees,
and clusters them based on a matrix of between-tree distances. Finally, it calculates a single representative tree for each cluster.

You can read the paper [here](https://dx.doi.org/10.1093/molbev/msw038)

## Installation

#### Preparing dependencies

If your system already has python 2.7, cython, numpy and a C++11-capable compiler (e.g. gcc >= 4.7), then you're ready to install.

The remaining python dependencies will be automatically installed during the build process.

#### External dependencies

To be able to build trees, treeCl needs to call on some external software. The choices are RAxML, PhyML, FastTree or PLL (using [pllpy](https://github.com/kgori/pllpy)). If any of these are installed, available in your path, and keep the standard names they were installed with, they should work.

#### Installing ``treeCl``
All remaining dependencies will be installed automatically using pip

    pip install treeCl



## Example Analysis
``` python

import treeCl

"""
The first point of call is the treeCl.Collection class. 
This handles loading your data, and calculating the trees 
and distances that will be used later.

This is how to load your data. This should be a directory
full of sequence alignments in fasta '*.fas' or phylip
'*.phy' formats. These can also be zipped using gzip or 
bzip2, treeCl will load them directly.
"""
c = treeCl.Collection(input_dir='input_dir', file_format='phylip')

"""
Now it's time to calculate some trees. The simplest way to 
do this is
"""
c.calc_trees()

"""
This uses RAxML to infer a tree for each alignment. We can 
pass arguments to RAxML using keywords.
"""
c.calc_trees(executable='raxmlHPC-PTHREADS-AVX',  # specify raxml binary to use
             threads=8,  # use multithreaded raxml
             model='PROTGAMMAWAGX',  # this model of evolution
             fast_tree=True)  # use raxml's experimental fast tree search option

"""
We can use PhyML instead of RAxML. Switching programs is 
done using a TaskInterface
"""

phyml = treeCl.tasks.PhymlTaskInterface()
c.calc_trees(task_interface=phyml)

"""
PhyML doesn't support multithreading, but treeCl can run 
multiple instances using JobHandlers
"""

threadpool = treeCl.parutils.ThreadpoolJobHandler(8)  # external software can be run in parallel
                                              # using a threadpool.
                                              
c.calc_trees(jobhandler=threadpool, task_interface=phyml)

"""
Trees are expensive to calculate. Results can be cached to disk, 
and reloaded.
"""
c.write_parameters('cache')
c = treeCl.Collection(input_dir='input_dir', param_dir='cache')

"""
Once trees have been calculated, we can measure all the 
distances between them. treeCl implements Robinson-Foulds (rf), 
weighted Robinson-Foulds (wrf), Euclidean (euc), and 
geodesic (geo) distances.
"""
dm = c.get_inter_tree_distances('geo')  

# Alternatively
processes = treeCl.parutils.ProcesspoolJobHandler(8)  # with pure python code, it is better to use processpools to parallelise for speed
dm = c.get_inter_tree_distances('geo', 
                                jobhandler=processes, 
                                batchsize=100)  # jobs are done in batches to
                                                # reduce overhead

"""
Hierarchical Clustering
"""
hclust = treeCl.Hierarchical(dm)
partition = hclust.cluster(3)  # partition into 3 clusters

# To use different linkage methods
from treeCl.clustering import linkage
partition = hclust.cluster(3, linkage.AVERAGE)
partition = hclust.cluster(3, linkage.CENTROID)
partition = hclust.cluster(3, linkage.COMPLETE)
partition = hclust.cluster(3, linkage.MEDIAN)
partition = hclust.cluster(3, linkage.SINGLE)
partition = hclust.cluster(3, linkage.WARD)  # default, Ward's method
partition = hclust.cluster(3, linkage.WEIGHTED)

"""
Spectral Clustering
"""
spclust = treeCl.Spectral(dm)
partition = spclust.cluster(3)

# Alternative calls
from treeCl.clustering import spectral, methods
spclust.cluster(3, algo=spectral.SPECTRAL, method=methods.KMEANS) # these are the defaults
spclust.cluster(3, algo=spectral.KPCA, method=methods.GMM) # alternatives use kernel PCA and a Gaussian Mixture Model

# Getting transformed coordinates
spclust.spectral_embedding(2) # spectral embedding in 2 dimensions
spclust.kpca_embedding(3) # kernel PCA embedding in 3 dimensions

"""
Multidimensional scaling
"""
mdsclust = treeCl.MultidimensionalScaling(dm)
partition = mdsclust.cluster(3)

# Alternatives: classical or metric MDS
from treeCl.clustering import mds
partition = mdsclust.cluster(3, algo=mds.CLASSICAL, method=methods.KMEANS)
partition = mdsclust.cluster(3, algo=mds.METRIC, method=methods.GMM)

# Getting transformed coordinates
mdsclust.dm.embedding(3, 'cmds')  # Classical MDS, 3 dimensions
mdsclust.dm.embedding(2, 'mmds')  # Metric MDS, 2 dimensions

"""
Score the result via likelihood
"""
raxml = treeCl.tasks.RaxmlTaskInterface()
sc = treeCl.Scorer(c, cache_dir='scorer', task_interface=raxml) 
sc.write_partition(partition)
results = sc.analyse_cache_dir(executable='raxmlHPC-PTHREADS-AVX', threads=8)

"""
Get the results
"""
# Get concatenated sequence alignments for each group
concats = [c.concatenate(grp) for grp in partition.get_membership()]
alignments = [conc.alignment for conc in concats]

# Get a list of the loci in each group
loci = sc.get_partition_members(partition)

# Get trees for each group
trees = sc.get_partition_trees(partition)

# Get full model parameters for each group
full_results = sc.get_partition_results(partition)  # same as returned by analyse_cache_dir

