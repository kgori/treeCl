#!/usr/bin/env python

"""
Tree DP means
"""

import numpy as np
import sys
from collections import defaultdict
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Partition
# from treeCl.lib.remote.datastructs import TrClTree

# class ExcessiveIndexing(Exception):

#     def __init__(self, index, length):
#         self.index = index
#         self.length = length

#     def __str__(self):
#         msg = ''

class DP_trees(object):

    """ Placeholder docstring """

    sumD_old = np.inf
    sumD_new = np.inf
    eps = sys.float_info.epsilon

    def __init__(self, scorer, lmbda, maxiter=1000):
        self.scorer = scorer
        self.lmbda = lmbda
        self.n = len(scorer.records)
        self.assignments = np.zeros(self.n, dtype=np.int)
        self.cluster_trees = []
        self.update_summary_trees()
        self.maxiter = maxiter


    @property
    def records(self):
        return self.scorer.records

    def distance(self, datapoint, cluster):
        return datapoint.tree.eucdist(cluster)

    def get_index_list(self, k):
        return np.where(self.assignments==k)[0]

    def get_summary_tree(self, k):
        ix = self.get_index_list(k)
        if not ix.any(): return None
        tree = self.scorer.add(tuple(ix))
        return tree

    def update_summary_trees(self):
        for i in range(self.assignments.max() + 1):
            summary_tree = self.get_summary_tree(i)
            self.set_ith_cluster_tree(i, summary_tree)

    def create(self, i):
        """
        Creates a new cluster and assigns current data to it
        """
        cluster = self.assignments.max() + 1
        tree = self.records[i].tree.copy()
        self.assignments[i] = cluster
        self.set_ith_cluster_tree(cluster, tree)

    def destroy(self, k):
        """
        Destroys a cluster, leaving members unassigned
        """
        ix = self.get_index_list(k)
        self.assignments[ix] = -1

    def assign_orphans(self):
        orphan_ix = self.get_index_list(-1)
        pass


    def set_ith_cluster_tree(self, i, tree):
        if len(self.cluster_trees) < (i+1):
            self.cluster_trees.append(None)
        if tree:
            self.cluster_trees[i] = tree.copy()

    def assignments_to_Partition(self):
        return Partition(self.assignments)

    def score(self):
        return self.scorer.score(self.assignments_to_Partition())

    def update_assignments(self): 
        for i, rec in enumerate(self.records):
            dists = [self.distance(rec, t) for t in self.cluster_trees]
            print dists
            k, min_dist = min(enumerate(dists), key=lambda x: x[1])

            if min_dist < self.lmbda:
                self.assignments[i] = k

            else: 
                self.create(i)


    def converged(self):
        """
        RANDOM PREDICATE PLACEHOLDER
        """
        if np.random.random() < 0.001:
            return True
        return False


    def run(self):
        itr = 0

        while not self.converged() and itr < self.maxiter:
            print 'iter',itr
            self.update_assignments()
            self.update_summary_trees()
            itr += 1

        return self.assignments

