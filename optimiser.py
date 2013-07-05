from treeCl import Collection, Scorer, Partition
from treeCl.lib.remote.utils import fileIO
from treeCl.lib.remote.externals.phyml import Phyml
from collections import defaultdict
import numpy as np
import random
import sys

EPS = 1e-8

class Optimiser(object):

    def __init__(self, nclusters, collection, datatype, tmpdir='/tmp', 
        initial_assignment=None):
        self.Collection = collection
        self.Scorer = Scorer(self.Collection.records, analysis='nj', 
                            datatype=datatype,
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


    def __reset_counts(self):
        self.done_worse = 0
        self.stayed_put = 0
        self.i = 0
        self.resets = 0

    def __status(self):
        return '{0} {1} {2}'.format(self.i, self.global_best_score,
            self.global_best_assignment)

    def initial_assignment(self, nclusters):
        return Partition(tuple(np.random.randint(nclusters, 
                size=len(self.Collection))))

    def get_cluster_trees(self, assignment):
        pvec = assignment.partition_vector
        index_dict = defaultdict(list)
        for (position, value) in enumerate(pvec):
            index_dict[value].append(position)
        tree_dict = {}
        for (k,v) in index_dict.items():
            tree_dict[k] = self.Scorer.concats[v]
        # trees = [self.Scorer.concats[i] for i in index_list]
        return tree_dict

    def move__(self, sample_size, assignment, nreassign=1, choose='max'):
        """
        MAKES A NEW PARTITION BY REASSIGNING RECORDS BETWEEN CLUSTERS
        How the ravel / unravel bit works:
        1:scores.ravel converts 2D scores array into 1D array
        2:argsort gives the indices that would sort the 1D array (ascending)
        3:np.unravel_index converts the 1D indices back into 2D indices,
          with 2D array size specified by scores.shape
        Result: list of n (i,j) pairs, where i is the sample index of the 
                test record (so Collection.records[sample[i]] is the record 
                itself), and j is the cluster it best fits.
        """
        sample = random.sample(range(len(self.Collection)), sample_size)
        cluster_trees = self.get_cluster_trees(assignment)
        scores = np.zeros((sample_size, self.nclusters))
        for i, record_index in enumerate(sample):
            rec = self.Collection.records[record_index]
            for j, tree in enumerate(cluster_trees):
                scores[i,j] = self.test(rec, tree)
        if choose=='max':
            moves = [np.unravel_index(i, scores.shape) 
                        for i in scores.ravel().argsort()[-nreassign:]]
        elif choose=='min':
            moves = [np.unravel_index(i, scores.shape) 
                        for i in scores.ravel().argsort()[:nreassign]]
        new_assignment = list(assignment.partition_vector)
        print moves
        for i,j in moves:
            new_assignment[sample[i]] = j+1 #because cluster number is in range
                                            # [1,x], and j is in range [0,x-1]
        return Partition(tuple(new_assignment)) 
        
    def move(self, sample_size, assignment, nreassign=1, choose='max'):
        """
        MAKES A NEW PARTITION BY REASSIGNING RECORDS BETWEEN CLUSTERS
        How the ravel / unravel bit works:
        1:scores.ravel converts 2D scores array into 1D array
        2:argsort gives the indices that would sort the 1D array (ascending)
        3:np.unravel_index converts the 1D indices back into 2D indices,
          with 2D array size specified by scores.shape
        Result: list of n (i,j) pairs, where i is the sample index of the 
                test record (so Collection.records[sample[i]] is the record 
                itself), and j is the cluster it best fits.
        """
        sample = random.sample(range(len(self.Collection)), sample_size)
        cluster_trees = self.get_cluster_trees(assignment)
        scores = np.zeros((sample_size, self.nclusters))
        for i, record_index in enumerate(sample):
            rec = self.Collection.records[record_index]
            for j, tree in enumerate(cluster_trees):
                scores[i,j] = self.test(rec, tree)
        
        new_clusters = scores.argmax(axis=1)
        M=scores/scores.sum(axis=1)[:, np.newaxis]
        reassignments = M.max(axis=1).argsort()[-nreassign:]

        new_assignment = list(assignment.partition_vector)

        for i in reassignments:
            new_assignment[sample[i]] = new_clusters[i]+1 #because cluster number is in range
                                            # [1,x], and new_clusters is in range [0,x-1]
        return Partition(tuple(new_assignment))     

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
        p.add_flag('-b', '0')   # no bootstraps
        p.add_flag('-m', model) # evolutionary model
        p.add_flag('-o', 'n')   # no optimisation 
        p.add_flag('-d', 'aa')  # datatype
        return p.run().score

    def optimise(self, 
        sample_size=10,
        nreassign=1,
        max_stayed_put=5, 
        max_resets=5, 
        max_done_worse=5, 
        max_iter=100):
        print 'Optimising - iter, best_score, partition'
        print self.__status()

        i = 0
        current_assignment = self.global_best_assignment
        while self.stayed_put < max_stayed_put:
            if self.resets == max_resets:
                print 'Reset limit reached ({0})'.format(max_resets)
                break
            if self.done_worse == max_done_worse:
                print 'wandered off, resetting...'
                self.resets += 1
                self.done_worse = 0
                current_assignment = self.global_best_assignment
            if i == max_iter:
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
            print self.__status()
            i += 1

        print self.__status()
        self.__reset_counts()


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(prog=fileIO.basename(__file__),
        description='Clustering optimiser')
    parser.add_argument('-n', '--nclusters')
    parser.add_argument('-f', '--format')
    parser.add_argument('-d', '--datatype')
    parser.add_argument('-p', '--filepath')
    parser.add_argument('-c', '--compression')
    parser.add_argument('-t', '--tmpdir')
