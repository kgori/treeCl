#!/usr/bin/env python

import os
import numpy as np
from collection import Scorer
from random import randint
from clustering import Partition
from distance_matrix import DistanceMatrix
from lib.remote.externals.phyml import Phyml


class EMTrees(object):
    def __init__(
        self,
        collection,
        nclusters,
        metric='euc',
        tmpdir=None,
    ):

        if not isinstance(nclusters, int) or nclusters <= 1:
            raise Exception('Need appropriate value for number of clusters.')

        self.nclusters = nclusters
        self.scorer = Scorer(collection.records, collection.analysis)  # Could check for entries
        self.datatype = collection.datatype
        self.metric = metric

        try:
            self.tmpdir
        except:
            self.tmpdir = collection.tmpdir

    def clusters_init(self):
        k = self.nclusters
        assignment = [0] * len(self.scorer.records)
        for i in range(k):
            assignment[np.random.randint(0, len(assignment))] = i + 1
        partition = Partition(assignment)
        clusters = [0] * k
        members = partition.get_membership()[1:]
        self.assign_clusters(clusters, members)
        for (index, record) in enumerate(self.scorer.records):
            scores = [self.ml(record, clusters[n]) for n in range(self.nclusters)]
            # print scores
            if assignment.count(assignment[index]) > 1 or assignment[index] == 0:
                assignment[index] = scores.index(max(scores)) + 1
        self.partition = Partition(assignment)
        self.L = self.scorer.score(self.partition)

    def random_partition(self):
        self.partition = Partition(tuple(np.random.randint(self.nclusters,
                                   size=len(self.scorer.records))))
        self.L = self.scorer.score(self.partition)

    def assign_clusters(self, clusters, members):
        for n in range(self.nclusters):
            if not clusters[n] or clusters[n].members != members[n]:
                clusters[n] = Cluster(members[n], self.scorer.records, self.scorer.analysis)

        return(clusters)

    def maximise(self, method):
        clusters = [0] * self.nclusters
        alg = getattr(self, method)
        count = 0

        while True:
            self.assign_clusters(clusters, self.partition.get_membership())
            assignment = list(self.partition.partition_vector)

            for (index, record) in enumerate(self.scorer.records):
                scores = [alg(record, clusters[n]) for n in range(self.nclusters)]
                # print scores
                if assignment.count(assignment[index]) > 1 or assignment[index] == 0:
                    assignment[index] = scores.index(max(scores)) + 1

            assignment = Partition(assignment)
            score = self.scorer.score(assignment)

            if score > self.L:
                self.L = score
                self.partition = assignment

            else:
                count += 1
                if count > 1: break  # Algorithm is deterministic so no need for more iterations

    def maximise_random(self, method):
        clusters = [0] * self.nclusters
        alg = getattr(self, method)
        count = 0
        sampled = []

        while True:
            self.assign_clusters(clusters, self.partition.get_membership())
            assignment = list(self.partition.partition_vector)

            index = randint(0, len(self.scorer.records) - 1)

            if index in sampled:
                continue
            else:
                record = self.scorer.records[index]
                sampled.append(index)

            scores = [alg(record, clusters[n]) for n in range(self.nclusters)]

            if assignment.count(assignment[index]) > 1 or assignment[index] == 0:
                assignment[index] = scores.index(max(scores)) + 1

            assignment = Partition(assignment)
            score = self.scorer.score(assignment)

            if score > self.L:
                self.L = score
                self.partition = assignment
                sampled = []
                count = 0
            else:
                count += 1
                if count == len(assignment): break

    def maximise_heuristic(self):
        clusters = [0] * self.nclusters
        sampled = []

        for i in range(1000):
            self.assign_clusters(clusters, self.partition.get_membership())
            assignment = list(self.partition.partition_vector)

            index = randint(0, len(self.scorer.records) - 1)

            record = self.scorer.records[index]
            sampled.append(index)

            lls = [self.ml(record, clusters[n]) for n in range(self.nclusters)]

            a = {'ll': max(lls)}
            a['n'] = lls.index(a['ll'])
            lls.pop(a['n'])

            b = {'ll': max(lls)}
            b['n'] = lls.index(b['ll'])

            a['p'] = np.maths.exp(a['ll'] - logsum(a['ll'], b['ll']))

            if np.random.uniform() > a['p']:
                choice = a['n']
            else:
                choice = b['n']

            if assignment.count(assignment[index]) > 1 or assignment[index] == 0:
                assignment[index] = choice + 1

            assignment = Partition(assignment)

            if i % 10 == 0:
                score = self.scorer.score(assignment)

                if score > self.L:
                    self.max_L = score
                    self.max_partition = assignment

    def dist(self, obj1, obj2):
        distance = DistanceMatrix([obj1.tree, obj2.tree], self.metric)[0][1]
        return(-distance)

    def ml(self, record, cluster, verbose=1):
        p = Phyml(record, tmpdir=self.tmpdir)
        input_tree = os.path.join(self.tmpdir, 'input_tree')
        cluster.tree.write_to_file(input_tree)
        p.add_tempfile(input_tree)
        p.add_flag('--inputtree', input_tree)
        p.add_flag('-o', 'r')  # Optimise only on substitutions`
        p.add_flag('-a', 'e')
        p.add_flag('-b', 0)
        p.add_flag('-c', 4)
        p.add_flag('--quiet', '')

        if self.datatype == 'protein':
            p.add_flag('-d', 'aa')
        elif self.datatype == 'dna':
            p.add_flag('-d', 'nt')

        score = p.run(verbosity=verbose).score
        return(score)


class Cluster(object):
    def __init__(self, members, records, analysis):
        self.members = tuple(members)
        self.records = [records[i] for i in self.members]
        self.scorer = Scorer(records, analysis)
        self.tree = self.scorer.add(self.members)


def logsum(loga, logb):
    # loga should be the larger
    b_a = 10**(logb - loga)
    return(loga + np.log10(1 + b_a))