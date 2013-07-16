#!/usr/bin/env python
from treeCl import Collection, Scorer, Partition
from treeCl.lib.remote.utils import fileIO
from treeCl.lib.remote.externals.phyml import Phyml
from collections import defaultdict
from analysis import Result
import operator
import numpy as np
import random
# import sys

EPS = 1e-8


class Optimiser(object):

    def __init__(self, nclusters, collection, tmpdir='/tmp',
                 initial_assignment=None):
        self.Collection = collection

        if not self.Collection.records[0].tree:
            print 'Calculating NJ trees for collection...'
            self.Collection.calc_NJ_trees()

        self.datatype = collection.datatype
        self.Scorer = Scorer(self.Collection.records, analysis='nj',
                             datatype=self.datatype,
                             tmpdir=tmpdir)

        if initial_assignment is None:
            initial_assignment = self.initial_assignment(nclusters)

        self.nclusters = nclusters
        self.tmpdir = tmpdir
        print 'Calculating initial scores...'
        self.global_best_score = self.Scorer.score(initial_assignment)
        self.global_best_assignment = initial_assignment
        self.done_worse = 0
        self.stayed_put = 0
        self.i = 0
        self.resets = 0

    def _reset_counts(self):
        self.done_worse = 0
        self.stayed_put = 0
        self.i = 0
        self.resets = 0

    def _status(self):
        return '{0} {1} {2}'.format(self.i, self.global_best_score,
                                    self.global_best_assignment)

    def initial_assignment(self, nclusters):
        return Partition(tuple(np.random.randint(nclusters,
                         size=len(self.Collection))))

    def get_clusters(self, assignment=None):
        p = assignment or self.global_best_assignment
        pvec = p.partition_vector
        index_dict = defaultdict(list)
        for (position, value) in enumerate(pvec):
            index_dict[value].append(position)
        return index_dict

    def get_cluster_trees(self, assignment, index_dict=None):
        index_dict = (index_dict or self.get_clusters(assignment))
        tree_dict = {}
        for (k, v) in index_dict.items():
            if not tuple(v) in self.Scorer.concats:
                self.Scorer.add(tuple(v))
            tree_dict[k] = self.Scorer.concats[tuple(v)]
        return tree_dict

    def score_sample(self, sample, assignment):
        """
        !! changed to simply SCORE a PRE-MADE SAMPLE
        sample_size:int, assignment:Partition object
        Calculates score m*n score matrix, where m is number of alignments
        in the sample, and n in the number of clusters encoded in the
        assignment (==Partition object)
        """
        # sample = random.sample(range(len(self.Collection)), sample_size)
        cluster_trees = self.get_cluster_trees(assignment)
        scores = np.zeros((len(sample), self.nclusters))
        for i, record_index in enumerate(sample):
            rec = self.Collection.records[record_index]
            for j, tree in cluster_trees.items():
                scores[i, j-1] = self.test(rec, tree)
        return (scores)

    def make_new_assignment(self, sample, scores, assignment, nreassign=1, choose='max'):
        """
        MAKES A NEW PARTITION BY REASSIGNING RECORDS BETWEEN CLUSTERS
        """

        new_clusters = scores.argmax(axis=1)
        M = scores/scores.sum(axis=1)[:, np.newaxis]
        if choose == 'max':
            reassignments = M.max(axis=1).argsort()[-nreassign:]
        elif choose == 'min':
            reassignments = M.min(axis=1).argsort()[:nreassign]

        new_assignment = list(assignment.partition_vector)

        for i in reassignments:
            new_assignment[sample[i]] = new_clusters[i]+1  # because cluster number is in range
                                            # [1,x], and new_clusters is in range [0,x-1]
        return Partition(tuple(new_assignment))

    def move(self, sample_size, assignment, nreassign=1, choose='max'):
        """
        !! now generates own sample and passes to scores
        wraps self.score_sample + self.new_assignment
        """
        sample = random.sample(range(len(self.Collection)), sample_size)
        scores = self.score_sample(sample, assignment)
        return self.make_new_assignment(sample, scores, assignment, nreassign,
                                        choose)

    def merge(self, assignment, label1, label2):
        pvec = ((x if x != label1 else label2)
                for x in assignment.partition_vector)
        return Partition(tuple(pvec))

    def merge_closest(self, assignment):
        print 'Finding clusters to merge...'
        clusters = self.get_clusters(assignment)
        best_score = -np.inf

        for i in clusters:
            for j in clusters:
                if i == j:
                    continue
                test_assignment = self.merge(assignment, i, j)
                score = self.Scorer.score(assignment)
                if score > best_score:
                    best_score = score
                    best_assignment = test_assignment

        return(best_assignment)

    def split(self, k, assignment=None):
        """
        Function to split cluster based on least representative alignment
        """
        assignment = assignment or self.global_best_assignment
        members = self.get_clusters(assignment)[k]
        tree_scores = {}
        cluster_record = self.Scorer.concatenate(members)
        print 'Calculating alignment scores...'

        for i in members:
            tree = self.Collection.records[i].tree
            tree_scores[i] = self.test(cluster_record, tree)

        print tree_scores

        seed, min_score = min(tree_scores.iteritems(), key=operator.itemgetter(1))
        print 'Splitting on {0}.'.format(seed)

        new_assignment = list(assignment.partition_vector)
        self.nclusters += 1
        new_assignment[seed] = self.nclusters
        print 'New Partition: {0}'.format(new_assignment)
        print 'Assigning to new partition...'

        new_assignment = Partition(new_assignment)
        scores = self.score_sample(members, new_assignment)
        assignment = self.make_new_assignment(members, scores, new_assignment, nreassign=len(members), choose='max')
        print 'Returning: {0}'.format(assignment)

        return assignment

    def split_max_var(self, assignment):
        clusters = self.get_clusters(assignment)
        var_dict = {}
        for k in clusters.keys():
            var_dict[k] = self.var(clusters[k])

        cluster_to_split, var = max(clusters.iteritems(), key=operator.itemgetter(1))

        return self.split(cluster_to_split, assignment)

    def test(self, record, tree, model='WAG'):
        """
        TESTS AN ALIGNMENT AGAINST A TREE TOPOLOGY
        """
        alignment_file = record.write_phylip('{0}/tmp_alignment.phy'.format(
            self.tmpdir), interleaved=True)
        newick_file = tree.write_to_file('{0}/tmp_tree.nwk'.format(self.tmpdir))
        p = Phyml(record)
        p.add_tempfile(alignment_file)
        p.add_tempfile(newick_file)
        p.add_flag('-i', alignment_file)
        p.add_flag('-u', newick_file)
        p.add_flag('-b', '0')    # no bootstraps
        p.add_flag('-m', model)  # evolutionary model
        p.add_flag('-o', 'n')    # no optimisation
        p.add_flag('-d', 'aa')   # datatype
        return p.run().score

    def var(self, members):
        score = self.Scorer.add(tuple(members)).score
        records = [self.Collection.records[i] for i in members]
        total_length = sum([r.seqlength for r in records])
        return(score / total_length)

    def optimise(self,
                 sample_size=10,
                 nreassign=1,
                 max_stayed_put=5,
                 max_resets=10,
                 max_done_worse=5,
                 max_iter=100):
        print 'Optimising - iter, best_score, partition'
        print self._status()

        current_assignment = self.global_best_assignment
        while True:
            if self.stayed_put > max_stayed_put:
                print 'stayed put too many times ({0})'.format(max_stayed_put)
                break
            if self.resets == max_resets:
                print 'Reset limit reached ({0})'.format(max_resets)
                break
            if self.done_worse == max_done_worse:
                print 'wandered off, resetting...'
                self.resets += 1
                self.done_worse = 0
                current_assignment = self.global_best_assignment
            if self.i == max_iter:
                print 'max iterations reached'
                break

            new_assignment = self.move(sample_size, current_assignment,
                                       nreassign)
            score = self.Scorer.score(new_assignment)
            print score, new_assignment

            if score > self.global_best_score:
                self.global_best_score = score
                self.global_best_assignment = new_assignment
                self.stayed_put = 0
                self.done_worse = 0
            elif np.abs(score - self.global_best_score) < EPS:
                self.stayed_put += 1
                self.done_worse = 0
            else:
                self.stayed_put = 0
                self.done_worse += 1
            print self._status()
            self.i += 1

        print self._status()
        self._reset_counts()


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(prog=fileIO.basename(__file__),
                                     description='Clustering optimiser')
    parser.add_argument('-n', '--nclusters', type=int)
    parser.add_argument('-f', '--format', default='phylip')
    parser.add_argument('-d', '--datatype', default='protein')
    parser.add_argument('-p', '--filepath', default='./')
    parser.add_argument('-c', '--compression', default=None)
    parser.add_argument('-t', '--tmpdir', default='/tmp/')
    # Collect all args for optimse and parse them later?
    parser.add_argument('-r', '--nreassign', default=10, type=int)
    parser.add_argument('-s', '--sample', default=10, type=int)

    args = parser.parse_args()

    print args

    c = Collection(input_dir=args.filepath,
                   compression=args.compression,
                   file_format=args.format,
                   datatype=args.datatype)

    o = Optimiser(args.nclusters, c)
    o.optimise(nreassign=args.nreassign, sample_size=args.sample)
    r = Result(o.global_best_score, o.global_best_assignment, o.Scorer.history)
    r.print_table()


"""
Changelog

renamed name-mangled classes
de-linted
inherited datatype from collection
"""
