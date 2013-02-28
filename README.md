# ``TreeCl`` - Phylogenetic Tree Clustering

TreeCl is a python package designed to take a collection 
of sequence alignments (of genes), infer their phylogenetic trees, and
perform clustering based on a distance matrix of their tree 
topologies.

The purpose of this is to establish whether there is any underlying structure
to the data.

### Dependencies:

- Python:
    - numpy
    - scipy
    - dendropy
    - scikit-learn (sklearn)

- External:
    - ALF  
    - Darwin
    - PhyML
    - TreeCollection
    - [GTP][]

[GTP]: <http://dx.doi.org/10.1109/TCBB.2010.3> "Owen, Megan, and J Scott Provan. 2011"