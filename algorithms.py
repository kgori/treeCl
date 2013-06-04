#!/usr/bin/env python

from collection import Collection, Scorer
from random import randint
from clustering import Partition
from distance_matrix import DistanceMatrix
from lib.remote.externals.phyml import Phyml

class emtrees(object): 

    # Not sure if scorer required
    # Possible end criteria?

    def __init__(
        self, 
        collection, 
        nclusters, 
        metric = 'euc', 
        ): 

        if not isinstance(nclusters, int) or nclusters <= 1:
            raise Exception('Need appropriate value for number of clusters.')

        self.nclusters = nclusters
        self.scorer = Scorer(collection.records, collection.analysis) # Could check for entries
        self.datatype = collection.datatype
        self.tmpdir = collection.tmpdir
        self.metric = metric

    def assign_partition(self):
        # Collapse equivalent trees?
        k = self.nclusters
        self.partition = Partition( [randint(1,k) for rec in self.scorer.records] )
        self.L = self.scorer.score(self.partition)

    def assign_clusters(self,clusters,partition):
        for n in range(self.nclusters):
            members = partition.get_membership()[n]
            if not clusters[n] or clusters[n].members != members:
                clusters[n] = Cluster(members, self.scorer.records, self.scorer.analysis)

        return(clusters)

    def maximise(self, method): # Needs to change
        self.assign_partition()
        clusters = [0] * self.nclusters
        alg = getattr(self,method)
        count = 0
        # clusters = [ Cluster(self.partition.get_membership()[n], self.scorer,records, self.scorer.analysis) for n in range(self.clusters)]

        while True:
            self.assign_clusters(clusters,self.partition)
            assignment = list(self.partition.partition_vector)

            for (index, record) in enumerate(self.scorer.records):
                scores = [ alg(record, clusters[n]) for n in range(self.nclusters) ]
                print scores
                if assignment.count(assignment[index]) > 1:
                    assignment[index] = scores.index(max(scores)) + 1

            assignment = Partition(assignment)
            score = self.scorer.score(assignment)

            if score > self.L:
                self.L = score
                self.partition = assignment

            else: 
                count += 1
                if count > 1: break # Algorithm is deterministic so no need for more iterations

    def dist(self, obj1, obj2):
        distance = DistanceMatrix( [obj1.tree, obj2.tree], self.metric)[0][1] 
        return(-distance)

    def ml(self, record, cluster, verbose=1):
        p = Phyml(record, tmpdir=self.tmpdir)
        cluster.tree.write_to_file('test_tree')
        p.add_tempfile('test_tree')
        p.add_flag('--inputtree', 'test_tree') # Need tempdir????
        p.add_flag('-o', 'r') # Optimise only on substitutions`
        p.add_flag('-a', 'e')
        p.add_flag('-b', 0)
        p.add_flag('-c', 4)
        if self.datatype == 'protein':
            p.add_flag('-d', 'aa') # set data type - could inherit from Collection object?
        elif self.datatype == 'dna':
            p.add_flag('-d', 'nt')
        tree = p.run(verbosity=verbose)
        return(tree.score)

class Cluster(object):
    def __init__(self, members, records, analysis):
        self.members = tuple(members)
        self.records = [ records[i] for i in self.members ]
        self.scorer = Scorer(records, analysis)
        self.tree = self.scorer.add(self.members)
