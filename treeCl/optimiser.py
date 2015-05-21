#!/usr/bin/env python
from __future__ import print_function

# standard library
from collections import defaultdict
import copy
import operator
import os
import random

# third party
import numpy as np

# treeCl
from collection import Collection, Scorer
from treeCl import tasks
from treeCl.parutils import sequential_map
from partition import Partition
from constants import EPS, NEGINF, TMPDIR, ANALYSES
from errors import optioncheck
from utils import fileIO, pll_helpers


def get_clusters(assignment):
    pvec = assignment.partition_vector
    index_dict = defaultdict(list)
    for (position, value) in enumerate(pvec):
        index_dict[value].append(position)
    return index_dict

def tr2dst(nwk):
    t = treeCl.Tree(nwk)
    taxa = t.patristic._pat_dists.keys()
    dm = np.zeros((len(taxa), len(taxa)))
    for i in range(len(taxa)):
        for j in range(i+1, len(taxa)):
            try:
                dst = t.patristic._pat_dists[taxa[i]][taxa[j]]
            except KeyError:
                dst = t.patristic._pat_dists[taxa[j]][taxa[i]]
            dm[i,j]=dm[j,i]=dst
    return treeCl.DistanceMatrix(dm, names=[taxon.label for taxon in taxa])

def read_model(al):
    model = al.parameters.partitions.model
    if model is None:
        if al.is_dna():
            model = 'DNA'
        else:
            model = 'LG'
    return model

def get_distance_based_likelihood(al, nwk):
    seqnames = al.get_names()
    al_dists = treeCl.DistanceMatrix(np.array(al.parameters.partitions.distances), names=seqnames)
    al_vars = treeCl.DistanceMatrix(np.array(al.parameters.partitions.variances), names=seqnames)
    tree_dists = tr2dst(nwk)
    lnl = 0
    for (seq1, seq2) in itertools.combinations(seqnames, 2):
        err = tree_dists.df[seq1][seq2] - al_dists.df[seq1][seq2]
        sdev = np.sqrt(al_vars.df[seq1][seq2])
        lnl += ss.norm.logpdf(err, scale=sdev)
    return lnl

def get_markov_based_likelihood(al, nwk):
    model = al.parameters.partitions.model
    if model is None and al.is_protein():
        model = 'LG08'
    al.set_substitution_model(model)
    freqs = al.parameters.partitions.frequencies
    al.set_frequencies(freqs)
    alpha = al.parameters.partitions.alpha
    if alpha:
        al.set_gamma_rate_model(4, alpha)
    else:
        al.set_constant_rate_model()
    if al.is_dna():
        rates = al.parameters.partitions.rates
        al.set_rates(rates, 'acgt')
    al.initialise_likelihood(nwk)
    lk = al.get_likelihood()
    # Hack to clear likelihood and free memory
    al.set_substitution_model(model)
    return lk

def get_markov_based_likelihood2(al, nwk):
    alignment_file, to_delete = al.get_alignment_file(as_phylip=True)
    if to_delete:
        filelist = [alignment_file]
    else:
        filelist = []
    model = read_model(al)
    partitions = '{}, al = 1 - {}'.format(model, len(al))
    freqs = al.parameters.partitions.frequencies
    alpha = al.parameters.partitions.alpha
    rates = al.parameters.partitions.rates if al.is_dna() else None
    with treeCl.utils.fileIO.TempFileList(filelist):
        inst = pllpy.pll(alignment_file, qfile, nwk, 1, 12345)
        inst.set_frequencies(freqs, 0, False)
        inst.set_alpha(alpha, 0, False)
        if rates: inst.set_rates(rates, 0, False)
    return inst.get_likelihood()

def get_markov_based_likelihood3(al, treelist):
    alignment_file, to_delete = al.get_alignment_file(as_phylip=True)
    if to_delete:
        filelist = [alignment_file]
    else:
        filelist = []
    model = al.parameters.partitions.model
    if model is None and al.is_protein():
        model = 'LG'
    partitions = '{}, al = 1 - {}'.format(model, len(al))
    with treeCl.utils.fileIO.TempFileList(filelist):
        inst = pllpy.pll(alignment_file, partitions, treelist[0], 1, 12345)
        freqs = al.parameters.partitions.frequencies
        alpha = al.parameters.partitions.alpha
        inst.set_frequencies(freqs, 0, False)
        inst.set_alpha(1, 0, False)
        results = [inst.get_likelihood()]
        for nwk in treelist[1:]:
            inst.set_tree(nwk)
            results.append(inst.get_likelihood())
    return results

def load_instances(collection):
    lst = []
    for al in collection:
        alignment_file, to_delete = al.get_alignment_file(as_phylip=True)
        if to_delete:
            filelist = [alignment_file]
        else:
            filelist = []
        with fileIO.TempFileList(filelist):
            model = read_model(al)
            partitions = '{}, {} = 1 - {}'.format(model, al.name, len(al))
            freqs = al.parameters.partitions.frequencies
            alpha = al.parameters.partitions.alpha
            rates = al.parameters.partitions.rates if al.is_dna() else None
            ml_tree = al.parameters.ml_tree
            inst = pll_helpers.create_instance(alignment_file, partitions, ml_tree, 1)
            inst.set_frequencies(freqs, 0, False)
            inst.set_alpha(alpha, 0, False)
            if rates: inst.set_rates(rates, 0, False)
            lst.append(inst)
    return lst

def _check_for_required_info(collection):
    """
    Raises an exception if required info is missing for an alignment.
    Info required is:
      ml tree (al.parameters.ml_tree)
      ml parameters (al.parameters.partitions.{frequencies,rates,alpha}
      likelihood model (al.parameters.partitions.model)
    :param collection:
    :return:
    """
    for al in collection:
        if not hasattr(al, 'parameters'):
            raise AttributeError('No parameters for alignment {}'.format(al.name))
        if not hasattr(al.parameters, 'ml_tree'):
            raise AttributeError('No ML tree for alignment {}'.format(al.name))
        if al.parameters.ml_tree is None:
            raise AttributeError('No ML tree for alignment {}'.format(al.name))
        if not hasattr(al.parameters, 'partitions'):
            raise AttributeError('No partition information for alignment {}'.format(al.name))
        if not al.parameters.partitions:
            raise AttributeError('No partition information for alignment {}'.format(al.name))



class EM(object):

    def __init__(self, collection):
        """
        Sets data
        """
        _check_for_required_info(collection)
        self.collection = collection

    
    def get_likelihood(self, i, instances, tree):
        """
        Return the log likelihood == logP(alignment[i] | tree, parameters)
        :param i: index of the alignment under test
        :param tree: newick string of the tree under test
        :return: log-likelihood (float)
        """
        inst = instances[i]
        inst.set_tree(tree)
        return inst.get_likelihood()

    def get_avg_per_site_likelihood(self, i, instances, tree):
        return self.get_likelihood(i, instances, tree) / len(self.collection[i])

    def initialise_groups(self):
        pass

    def mstep(self):
        """
        In general, the process of maximising the likelihood by tuning the partition parameter
        Uses the following strategies:
          single reassignment
          batch reassignment
          split-and-merge
        :return:
        """
        pass

    def estep(self, partition):
        """
        In general, estimating the likelihood of the data using the current partition parameter value
        This involves computing new trees based on the partition assignments. Theoretically, this
        should be done by maximising the likelihood of each tree using the full set of parameters,
        but some shortcuts are taken:
          tree search: least-square trees + rough ML branch length optimisation
          model parameters: fixed at pre-estimated values

        :return:
        """
        memberships = partition.get_membership()
        args = []
        for mbl in memberships:
            args.append(self.collection.get_tree_collection_strings(mbl))
        msg = 'ESTEP'
        trees = [result['tree'] for result in sequential_map(tasks.minsq_task, args, msg)]


class Optimiser(object):
    def __init__(self, nclusters, collection, tmpdir=TMPDIR,
                 analysis='nj', initial_assignment=None, scorer=None):
        optioncheck(analysis, ANALYSES + ['tc', 'TreeCollection'])
        if self.analysis == 'tc':
            self.analysis = 'TreeCollection'
        else:
            self.analysis = analysis

        self.Collection = collection

        if not self.Collection.records[0].tree:
            print('Calculating {} trees for collection...'.format(analysis))
            self.Collection.calc_NJ_trees()

        self.datatype = collection.datatype
        if scorer is not None and isinstance(scorer, Scorer):
            self.scorer = scorer
        else:
            self.scorer = Scorer(self.Collection)

        self.nclusters = nclusters
        self.tmpdir = tmpdir

        print('Calculating initial scores...')
        if initial_assignment is None:
            initial_assignment = Partition(tuple([0] * len(collection)))
            # initial_assignment = self.random_partition(nclusters)

        self.global_best_scores = {}
        self.global_best_assignments = {}
        self.global_best_scores[self.nclusters] = self.scorer.score(
            initial_assignment, history=True)
        self.global_best_assignments[self.nclusters] = initial_assignment

        self.done_worse = 0
        self.stayed_put = 0
        self.i = 0
        self.resets = 0
        self.merges = 0

    def _reset_counts(self):
        self.done_worse = 0
        self.stayed_put = 0
        self.i = 0
        self.resets = 0

    def status(self, current_assignment, details=None):
        iter_ = self.i
        n = len(current_assignment)
        curr_score = self.scorer.score(current_assignment, history=False)
        best_score = self.global_best_scores[n]
        details = ('\t' + str(details) if details is not None else '')

        return 'Iter:{0}\tNclusters:{1}\tCurrent\tscore:{2}\tBest score:{3}{4}'.format(
            iter_, n, curr_score, best_score, details)

    def random_partition(self, nclusters):
        return Partition(tuple(np.random.randint(nclusters,
                                                 size=len(self.Collection))))

    def update(self, assignment):
        """
        method for working interactively and keeping nclusters correct
        """
        nclusters = len(assignment)  # len(assignment) == number of clusters
        best_score = self.global_best_scores.get(nclusters, NEGINF)
        curr_score = self.scorer.score(assignment, history=False)
        if (curr_score - best_score) > EPS:
            self.global_best_assignments[nclusters] = assignment
            self.global_best_scores[nclusters] = self.scorer.score(assignment,
                                                                   history=False)

    def get_cluster_trees(self, assignment, index_dict=None):
        index_dict = (index_dict or get_clusters(assignment))
        tree_dict = {}
        for (k, v) in index_dict.items():
            if not tuple(v) in self.scorer.concats:
                self.scorer.add(tuple(v))
            tree_dict[k] = self.scorer.concats[tuple(v)]
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
        scores = np.zeros((len(sample), len(cluster_trees)))
        for i, record_index in enumerate(sample):
            rec = self.Collection.records[record_index]
            for j, tree in cluster_trees.items():
                scores[i, j - 1] = self.test(rec, tree)
        return scores

    def constrain_assignment(self, assignment, nclusters=None):
        """
        Constrain the assignment to have self.nclusters clusters
        """

        if nclusters is None:
            nclusters = self.nclusters
        if (nclusters < 1) or (nclusters > len(self.Collection)):
            raise ValueError('Invalid number of clusters: {}'.format(nclusters))
        while len(assignment.get_membership()) > nclusters:
            assignment = self.merge_closest(assignment)
        while len(assignment.get_membership()) < nclusters:
            assignment = self.split_search(assignment)
        return assignment

    def make_new_assignment(self, sample, scores, assignment, nreassign=1,
                            choose='max'):
        """
        MAKES A NEW PARTITION BY REASSIGNING RECORDS BETWEEN CLUSTERS
        """
        optioncheck(choose, ('max', 'min'))
        new_clusters = scores.argmax(axis=1)
        M = scores / scores.sum(axis=1)[:, np.newaxis]
        if choose == 'max':
            reassignments = M.max(axis=1).argsort()[-nreassign:]
        else:
            reassignments = M.min(axis=1).argsort()[:nreassign]

        new_assignment = list(assignment.partition_vector)

        for i in reassignments:
            new_assignment[sample[i]] = new_clusters[i] + 1
            # because cluster number is in range
            # [1,x], and new_clusters is in range [0,x-1]

        return Partition(tuple(new_assignment))

    def move(self, sample_size, assignment, nreassign=1, choose='max',
             sampled=None):
        """
        !! now generates own sample and passes to scores
        wraps self.score_sample + self.new_assignment
        """

        if sampled is None:
            sampled = list()

        unsampled = set(range(len(self.Collection))) - set(sampled)

        if len(unsampled) > 0:
            if sample_size > len(unsampled):
                sample = unsampled
            else:
                sample = random.sample(unsampled, sample_size)

            self.sampled.extend(sample)
            scores = self.score_sample(sample, assignment)
            assignment = self.make_new_assignment(sample, scores, assignment,
                                                  nreassign, choose)
        return assignment

    def merge(self, assignment, label1, label2):
        pvec = ((x if x != label1 else label2)
                for x in assignment.partition_vector)
        return Partition(tuple(pvec))

    def merge_closest(self, assignment):
        print('Finding clusters to merge...')
        clusters = get_clusters(assignment)
        best_score = NEGINF
        merging = [None, None]

        for i in clusters:
            for j in clusters:
                # print('i = {}, j = {}'.format(i, j))
                if i >= j:
                    continue
                print('Testing Clusters {0} and {1}'.format(i, j))
                test_assignment = self.merge(assignment, i, j)
                self.update(test_assignment)
                score_value = self.scorer.score(test_assignment, history=False)

                if score_value > best_score:
                    merging[0] = i
                    merging[1] = j
                    best_score = score_value
                    best_assignment = test_assignment

        print('Merging clusters {0} and {1}'.format(*merging))
        print('Best assignment: {0}'.format(best_assignment))
        return best_assignment

    def split(self, k, assignment, verbosity=1):
        """
        Function to split cluster based on least representative alignment
        """
        if verbosity > 1:
            print(assignment)
        members = get_clusters(assignment)[k]
        if len(members) == 1:
            return assignment
        elif len(members) == 2:
            new_partition_vector = list(assignment.partition_vector)
            new_partition_vector[members[0]] = max(assignment.partition_vector) + 1
            new_assignment = Partition(new_partition_vector)
            return new_assignment

        tree = self.get_cluster_trees(assignment)[k]
        alignment_scores = {}
        if verbosity > 0:
            print('Calculating alignment scores...')

        for i in members:
            r = self.Collection.records[i]
            alignment_scores[i] = self.test(r, tree) / float(r.seqlength)
            # per-site likelihood

        seed, min_score = min(alignment_scores.iteritems(),
                              key=operator.itemgetter(1))
        print('Splitting on {0}.'.format(seed + 1))  # convert to 1-based indexing

        new_assignment = list(assignment.partition_vector)
        new_assignment[seed] = max(assignment.partition_vector) + 1
        if verbosity > 1:
            print('New Partition: {0}'.format(new_assignment))
        if verbosity > 0:
            print('Assigning to new partition...')

        new_assignment = Partition(new_assignment)
        scores = self.score_sample(members, new_assignment)
        assignment = self.make_new_assignment(members, scores, new_assignment,
                                              nreassign=len(members))
        if verbosity > 1:
            print('Returning: {0}'.format(assignment))

        return assignment

    def split_max_var(self, assignment):
        clusters = get_clusters(assignment)
        var_dict = {}

        for k in clusters.keys():
            var_dict[k] = self.var(clusters[k])

        print(var_dict)

        cluster_to_split, var = max(clusters.iteritems(),
                                    key=operator.itemgetter(1))

    def split_search(self, assignment, update=True):
        clusters = get_clusters(assignment)
        k = len(assignment)
        best_score = NEGINF

        for i in clusters:
            print('i: {0}'.format(i))
            test_assignment = self.split(i, assignment)
            # score = self.scorer.score(test_assignment)
            if len(test_assignment) == k + 1:
                curr_score = self.scorer.score(test_assignment, history=False)
                self.update(test_assignment)
            else:
                curr_score = -np.Inf
                print('Something has gone wrong')
            print(test_assignment)
            print(curr_score)

            if curr_score > best_score:
                best_score = curr_score
                best_assignment = test_assignment

        return best_assignment

    def test(self, record, tree, model=None):
        """
        TESTS AN ALIGNMENT AGAINST A TREE TOPOLOGY
        """
        tmp_record = copy.deepcopy(record)

        # if tree label set and record label set don't match
        header_set = set(tmp_record.headers)
        extra_in_tree = tree.labels - header_set
        extra_in_record = header_set - tree.labels

        if extra_in_tree:
            for lab in extra_in_tree:
                tmp_record.headers.append(lab)
                tmp_record.sequences.append(''.join(['-'] * tmp_record.seqlength))
            tmp_record.update()

        if extra_in_record:
            for lab in extra_in_record:
                i = tmp_record.headers.index(lab)
                tmp_record.headers = (tmp_record.headers[:i] +
                                      tmp_record.headers[i + 1:])
                tmp_record.sequences = (tmp_record.sequences[:i] +
                                        tmp_record.sequences[i + 1:])
            tmp_record.update()

        return tmp_alignment.likelihood(tree, self.tmpdir, fit_rates=True)
        # alignment_file = tmp_record.write_phylip('{0}/tmp_alignment.phy'.format(
        # self.tmpdir), interleaved=True)
        # newick_file = tree.write_to_file('{0}/tmp_tree.nwk'.format(self.tmpdir))
        # p = Phyml(tmp_record, self.tmpdir)
        # p.add_tempfile(alignment_file)
        # p.add_tempfile(newick_file)
        # p.add_flag('-i', alignment_file)
        # p.add_flag('-u', newick_file)
        # p.add_flag('-b', '0')    # no bootstraps
        # if tmp_record.datatype == 'dna':
        #     if model is None:
        #         model = 'GTR'
        #     p.add_flag('-m', model)
        #     p.add_flag('-d', 'nt')
        # else:
        #     if model is None:
        #         model = 'WAG'
        #     p.add_flag('-m', model)  # evolutionary model
        #     p.add_flag('-d', 'aa')   # datatype

        # p.add_flag('-o', 'n')    # no optimisation
        # return p.run().score

    def var(self, members):
        score = self.scorer.add(tuple(members)).score
        records = [self.Collection.records[i] for i in members]
        total_length = sum([r.seqlength for r in records])

        return score / total_length

    def optimise(self,
                 assignment,
                 nclusters=None,
                 update=True,
                 history=True,
                 sample_size=10,
                 nreassign=10,
                 max_stayed_put=25,
                 max_resets=5,
                 max_done_worse=5,
                 max_iter=1000):

        if nclusters is None:
            nclusters = self.nclusters

        assignment = self.constrain_assignment(assignment, nclusters)

        local_best_assignment = assignment
        local_best_score = self.scorer.score(local_best_assignment,
                                             history=False)
        current_assignment = local_best_assignment
        self.sampled = []

        print(self.status(current_assignment))

        while True:
            if self.stayed_put > max_stayed_put:
                print('stayed put too many times ({0})'.format(max_stayed_put))
                break
            if self.resets == max_resets:
                print('Reset limit reached ({0})'.format(max_resets))
                break
            if self.done_worse == max_done_worse:
                print('wandered off, resetting...')
                self.resets += 1
                self.done_worse = 0
                current_assignment = local_best_assignment
            if self.i == max_iter:
                print('max iterations reached')
                break

            new_assignment = self.move(sample_size, current_assignment,
                                       nreassign)
            new_assignment = self.constrain_assignment(new_assignment,
                                                       nclusters)
            score = self.scorer.score(new_assignment, history=history)
            self.update(new_assignment)

            if (score - local_best_score) > EPS:
                self.sampled = []
                local_best_score = score
                local_best_assignment = new_assignment
                self.stayed_put = 0
                self.done_worse = 0
                self.resets = 0
                print(self.status(new_assignment, '(Improved)'))
            elif np.abs(score - local_best_score) < EPS:
                self.stayed_put += 1
                self.done_worse = 0
                message = ('(No improvement - [{}/{}])'.format(self.stayed_put,
                                                               max_stayed_put))
                print(self.status(new_assignment, message))
            else:
                self.sampled = []
                # self.stayed_put = 0
                self.done_worse += 1
                message = '(Did worse - [{}/{}]'.format(self.done_worse,
                                                        max_done_worse)
                print(self.status(new_assignment, message))

            self.i += 1

        self._reset_counts()
        return local_best_assignment

    def optimise_with_variable_clusters(self,
                                        assignment,
                                        target_clusters,
                                        max_clusters,
                                        optimise_on_ascent=True,
                                        optimise_on_descent=True,
                                        update=True,
                                        **kwargs):

        if max_clusters < target_clusters:
            raise ValueError('max_clusters ({}) must be at least equal to '
                             'target_clusters ({})'.format(max_clusters, target_clusters))

        current_clusters = len(assignment)
        print('Optimising current assignment with {} clusters. Optimiser will '
              'ascend to {} clusters, and descend to a target of {} clusters'
              '.'.format(current_clusters, max_clusters, target_clusters))
        for n in range(current_clusters, max_clusters + 1):
            print("ASCENDING (optimisation:{}) -> Current target: "
                  "{} clusters".format(('ON' if optimise_on_ascent else 'OFF'),
                                       n))
            if optimise_on_ascent:
                assignment = self.optimise(assignment, nclusters=n, **kwargs)
            else:
                assignment = self.constrain_assignment(assignment, n)

        for n in range(max_clusters - 1, target_clusters - 1, -1):
            print('DESCENDING (optimisation:{}) -> Current target: {} '
                  'clusters'.format(('ON' if optimise_on_descent else 'OFF'),
                                    n))
            if optimise_on_descent:
                assignment = self.optimise(assignment, nclusters=n, **kwargs)
            else:
                assignment = self.constrain_assignment(assignment, n)

        return self.constrain_assignment(assignment, target_clusters)

    def write(self, filename):
        headers = ['Iteration', 'CPU Time', 'Likelihood', 'Partition',
                   'NClusters']
        output = [[i] + x + len(x[-1])
                  for (i, x) in enumerate(self.scorer.history)]

        with open(filename, 'w+') as file_:
            writer = csv.writer(file_, delimiter='\t', quoting=csv.QUOTE_NONE)
            writer.writerow(headers)
            writer.writerows(output)


def get_partition(clusters):
    seq = clusters if isinstance(clusters, dict) else range(len(clusters))
    length = sum([len(clusters[i]) for i in seq])
    pvec = [0] * length

    for k in seq:
        for i in clusters[k]:
            pvec[i] = k

    return Partition(tuple(pvec))


def get_partition_from_file(filename):
    with open(filename) as f:
        pvec = [int(x) for x in f.readline().split()]

    return Partition(tuple(pvec))


if __name__ == '__main__':
    import tempfile
    import string
    import csv
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Clustering optimiser')

    parser.add_argument('-n',
                        '--nclusters',
                        type=int)

    parser.add_argument('-f',
                        '--format',
                        default='phylip')

    parser.add_argument('-d',
                        '--datatype',
                        default='protein')

    parser.add_argument('-i',
                        '--input_dir',
                        default='./')

    parser.add_argument('-c',
                        '--compression',
                        default=None)

    parser.add_argument('-t',
                        '--tmpdir',
                        default='/tmp/')

    parser.add_argument('-r',
                        '--nreassign',
                        default=10, type=int)

    parser.add_argument('-w',
                        '--max_done_worse',
                        default=5)

    parser.add_argument('-s',
                        '--sample_size',
                        default=10,
                        type=int)

    parser.add_argument('-o',
                        '--output',
                        default=None)

    parser.add_argument('-m',
                        '--merge',
                        action='store_true',
                        help='Enable merge/splitting of clusters')

    parser.add_argument('-max',
                        type=int,
                        help='Maximum number of clusters to ascend to')

    parser.add_argument('-target',
                        type=int,
                        help='Target number of clusters to descend to')

    parser.add_argument('-p',
                        '--partition',
                        default=None,
                        help='Initial partitioning (e.g. from clustering), given as a file')

    parser.add_argument('--hierarchical',
                        default=None,
                        choices=['bottom_up', 'top_down'],
                        help=('use top-down/bottom-up hierarchical clustering to generate'
                              ' initial partitioning'))

    parser.add_argument('-q',
                        '--quit',
                        action='store_true',
                        help='Quit before optimising, just return initial partition and score')

    # ## GET COMMAND LINE ARGUMENTS
    args = parser.parse_args()

    if args.sample_size < args.nreassign:
        args.sample_size = args.nreassign
    opt_args = {
        'nreassign': args.nreassign, 'sample_size': args.sample_size,
        'max_done_worse': args.max_done_worse
    }

    new_tmpdir = tempfile.mkdtemp(prefix='tmpwrap_mgp_', dir=args.tmpdir)

    c = Collection(input_dir=args.input_dir,
                   compression=args.compression,
                   file_format=args.format)

    # Initial partitioning - either random or specified from a file
    if args.partition is not None:

        if not os.path.exists(args.partition):
            print('Partition file {0} not found'.format(args.partition))
            sys.exit(1)
        else:
            p = get_partition_from_file(args.partition)

        if not len(p) == len(c):
            print('Partition is of incorrect length '
                  '(expected {0}, got {1}'.format(len(c), len(p)))
            sys.exit(1)

        o = Optimiser(args.nclusters, c, tmpdir=new_tmpdir,
                      initial_assignment=p)

    else:
        o = Optimiser(args.nclusters, c, tmpdir=new_tmpdir)

    # Hierarchical clustering via likelihood
    if args.hierarchical is not None:
        if args.hierarchical == 'top_down':
            p = Partition(tuple([1] * len(c)))
        elif args.hierarchical == 'bottom_up':
            p = Partition(range(1, len(c) + 1))

        result = o.constrain_assignment(p, args.nclusters)
        # o.Scorer.clear_history()
        score = o.Scorer.score(result)
        o.global_best_assignments[args.nclusters] = result
        o.global_best_scores[args.nclusters] = score

    # Quit early
    if args.quit:
        pass

    else:
        if args.merge is True:
            o.optimise_with_merge(o.global_best_assignments[args.nclusters],
                                  update=True, **opt_args)
        else:
            o.optimise(o.global_best_assignments[args.nclusters],
                       update=True, **opt_args)

    def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    output_name = args.output or 'output_' + id_generator(6)
    o.write(output_name)

"""
Changelog

renamed name-mangled classes
de-linted
inherited datatype from collection
"""
