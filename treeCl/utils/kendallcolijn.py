from collections import namedtuple
from itertools import combinations
import numpy as np

class KendallColijn(object):
    """
    Data structure that stores info about a tree
    that is needed to compute the Kendall-Colijn
    tree  distance  metric  -  i.e. the  vectors 
    m and M
    """
    def __init__(self, tree):
        """
        Initialise the data structure, compute m and M.
        """
        info = self._precompute(tree)
        m, M = self._get_vectors(tree, info)
        self.little_m = m
        self.big_m = M

    def _precompute(self, tree):
        """
        Collect metric info in a single preorder traversal.
        """
        d = {}
        for n in tree.preorder_internal_node_iter():
            d[n] = namedtuple('NodeDist', ['dist_from_root', 'edges_from_root'])
            if n.parent_node:
                d[n].dist_from_root = d[n.parent_node].dist_from_root + n.edge_length
                d[n].edges_from_root = d[n.parent_node].edges_from_root + 1
            else:
                d[n].dist_from_root = 0.0
                d[n].edges_from_root = 0
        return d

    def _get_vectors(self, tree, precomputed_info):
        """
        Populate the vectors m and M.
        """
        little_m = []
        big_m = []

        leaf_nodes = sorted(tree.leaf_nodes(), key=lambda x: x.taxon.label)
        # inner nodes, sorted order
        for leaf_a, leaf_b in combinations(leaf_nodes, 2):
            mrca = tree.mrca(taxa=[leaf_a.taxon, leaf_b.taxon])
            little_m.append(precomputed_info[mrca].edges_from_root)
            big_m.append(precomputed_info[mrca].dist_from_root)
        
        # leaf nodes, sorted order
        for leaf in leaf_nodes:
            little_m.append(1)
            big_m.append(leaf.edge_length)

        return np.array(little_m), np.array(big_m)

    def get_vector(self, lbda=0.5):
        """
        The vector v is the weighted average of m and M.
        lbda, a.k.a. lambda, is the weighting parameter.
        """
        return (1-lbda)*self.little_m - lbda*self.big_m

    def get_distance(self, other, lbda=0.5):
        """
        Return the Euclidean distance between vectors v of
        two trees. Must have the same leaf set (too lazy to check).
        """
        return np.sqrt(((self.get_vector(lbda) - other.get_vector(lbda))**2).sum())
