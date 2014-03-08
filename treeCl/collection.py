#!/usr/bin/env python

# standard lib
import sys
import itertools
import timeit

# third party
from dendropy import TaxonSet

# treeCl
from datastructs.seq import concatenate
from datastructs.trcl_seq import TrClSeq
from datastructs.trcl_tree import TrClTree
from distance_matrix import DistanceMatrix
from software_interfaces.alf import simulate_from_record
from software_interfaces.DVscript import runDV
from software_interfaces.phyml import runPhyml, runLSFPhyml
from software_interfaces.treecollection import runTC
from utils import flatten_list
from errors import  OptionError, optioncheck, directorymake,\
    directorycheck
from utils import fileIO
from constants import TMPDIR, SORT_KEY, ANALYSES


class NoRecordsError(Exception):

    def __init__(self, file_format, input_dir, compression):
        self.file_format = file_format
        self.input_dir = input_dir
        self.compression = compression

    def __str__(self):
        msg = ('No records were found in {0} matching\n'
            '\tfile_format = {1}\n'
            '\tcompression = {2}'.format(self.input_dir,
                self.file_format, self.compression))
        return msg


class Collection(object):

    """ Call:

    c = Collection(inut_dir, file_format, datatype, tmpdir ...)
    c.calc_distances(), c.calc_TC_trees(), ...
    dm = c.distance_matrix('geo')
    cl = Clustering(dm)
    k = cl.spectral(4, prune='estimate', local_scale=7)
    p = Partition(k) """

    def __init__(
        self,
        records=None,
        input_dir=None,
        file_format='fasta',
        datatype=None,
        tmpdir=TMPDIR,
        calc_distances=False,
        compression=None,
        analysis=None,
        ):

        self.tmpdir = directorymake(tmpdir)

        if records:
            self.records = records
            self.datatype = datatype or records[0].datatype
            optioncheck(self.datatype, ['dna', 'protein'])
            for rec in records:
                rec.tmpdir = self.tmpdir

        elif input_dir:
            directorycheck(input_dir)
            self.datatype = optioncheck(datatype, ['dna', 'protein'])
            optioncheck(file_format, ['fasta', 'phylip'])
            self.records = self.read_files(input_dir, file_format, compression)

        else:
            print 'Provide a list of records, or the path to a set of alignments'

        if not self.records:
            raise NoRecordsError(file_format, input_dir, compression)

        if calc_distances:
            self.calc_distances()

        self.taxon_set = TaxonSet()

    def __len__(self):
        if getattr(self, 'records'):
            return len(self.records)
        return 0

    @property
    def records(self):
        return [self._records[i] for i in range(len(self._records))]

    @records.setter
    def records(self, records):

        for rec in records:
            rec.sanitise()
        self._records = dict(enumerate(records))

        # self.reverse_lookup = {v:k for (k,v) in self.records}

    @property
    def trees(self):
        return [self._records[i].tree for i in range(len(self._records))]

    def read_files(self, input_dir, file_format, compression=None):
        """ Get list of alignment files from an input directory *.fa, *.fas and
        *.phy files only

        Stores in self.files """

        optioncheck(compression, [None, 'gz', 'bz2'])

        if file_format == 'fasta':
            extensions = ['fa', 'fas', 'fasta']

        elif file_format == 'phylip':
            extensions = ['phy']

        if compression:
            extensions = ['.'.join([x, compression]) for x in extensions]

        files = fileIO.glob_by_extensions(input_dir, extensions)
        files.sort(key=SORT_KEY)

        return [TrClSeq(f, file_format=file_format, datatype=self.datatype,
                        name=fileIO.strip_extensions(f),
                        tmpdir=self.tmpdir)
                for f in files]

    def calc_distances(self, verbosity=0):
        for rec in self.records:
            runDV(rec, tmpdir=self.tmpdir, verbosity=verbosity)

    def calc_TC_trees(self, verbosity=0):
        self.analysis = 'TreeCollection'
        for rec in self.records:
            runTC(rec, self.tmpdir, verbosity=verbosity,
                taxon_set=self.taxon_set)
            rec.tree = TrClTree.cast(rec.tree)

    def calc_ML_trees(self, lsf=False, verbosity=0):
        """ Deprecated"""
        print 'deprecated: use calc_phyml_trees instead'
        self.analysis = 'ml'
        if lsf:
            trees = runLSFPhyml(self.records,
                                self.tmpdir,
                                analysis=self.analysis,
                                verbosity=verbosity)
            for rec, tree in zip(self.records, trees):
                rec.tree = TrClTree.cast(tree)
        else:
            for rec in self.records:
                runPhyml(rec, self.tmpdir, analysis=self.analysis,
                    verbosity=verbosity, taxon_set=self.taxon_set)
                rec.tree = TrClTree.cast(rec.tree)

    def calc_NJ_trees(self, lsf=False, analysis='nj', verbosity=0):
        print 'Deprecated:use calc_phyml_trees instead'
        self.calc_phyml_trees(self, lsf, analysis, verbosity)

    def calc_phyml_trees(self, lsf=False, analysis='nj', verbosity=0):
        optioncheck(analysis, ANALYSES)
        self.analysis = analysis
        if lsf:
            trees = runLSFPhyml(self.records,
                                self.tmpdir,
                                analysis=self.analysis,
                                verbosity=verbosity)
            for rec, tree in zip(self.records, trees):
                rec.tree = TrClTree.cast(tree)
        else:
            for rec in self.records:
                runPhyml(rec, self.tmpdir, analysis=analysis,
                         verbosity=verbosity, taxon_set=self.taxon_set)
                rec.tree = TrClTree.cast(rec.tree)
        if verbosity == 1:
            print

    def get_phyml_command_strings(self, analysis, tmpdir, verbosity=0):
        cmds = [runPhyml(rec, tmpdir, analysis=analysis,
                         verbosity=verbosity, taxon_set=self.taxon_set,
                         dry_run=True)
                for rec in self.records]
        return cmds

    def distance_matrix(self, metric, **kwargs):
        """ Generate a distance matrix from a fully-populated Collection """
        trees = [rec.tree for rec in self.records]
        return DistanceMatrix(trees, metric, tmpdir=self.tmpdir,
            **kwargs)

    def permuted_copy(self):
        lengths, names = zip(*[(rec.seqlength, rec.name) for rec in self.records])
        concat = concatenate(self.records)
        concat.shuffle()
        new_records = concat.split_by_lengths(lengths, names)
        return self.__class__(new_records)


class Scorer(object):

    """ Takes an index list, generates a concatenated SequenceRecord, calculates
    a tree and score """

    def __init__(
        self,
        records,
        analysis,
        lsf=False,
        max_guidetrees=10,
        tmpdir=None,
        datatype=None,
        verbosity=0,
        ):

        self.analysis = optioncheck(analysis, ANALYSES + ['TreeCollection'])
        self.max_guidetrees = max_guidetrees
        self.lsf = lsf
        self.records = records
        self.datatype = datatype or records[0].datatype
        self.verbosity = verbosity
        optioncheck(self.datatype, ['protein', 'dna'])
        self.tmpdir = tmpdir or records[0].tmpdir
        directorymake(self.tmpdir)
        self.concats = {}
        self.history = []
        self.populate_cache()

    def add_partition_list(self, partition_list):
        index_tuples = list(itertools.chain(*[partition.get_membership()
                                              for partition in partition_list]))
        missing = sorted(set(index_tuples).difference(self.concats.keys()))
        if self.lsf and not self.analysis == 'TreeCollection':
            supermatrices = [self.concatenate(index_tuple)
                             for index_tuple in index_tuples]
            trees = runLSFPhyml(supermatrices,
                                self.tmpdir,
                                analysis=self.analysis,
                                verbosity=self.verbosity)
            for tree in trees:
                tree = TrClTree.cast(tree)
            for index_tuple, tree in zip(missing, trees):
                self.concats[index_tuple] = tree
        else:
            for index_tuple in missing:
                self.add(index_tuple)

    def add(self, index_tuple):
        """ Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree """

        if index_tuple in self.concats:
            return self.concats[index_tuple]

        concat = self.concatenate(index_tuple)

        if self.analysis == 'TreeCollection':
            guidetrees = [self.records[n].tree for n in
                          index_tuple][:self.max_guidetrees]
            tree = TrClTree.cast(runTC(concat, self.tmpdir, guidetrees,
                                verbosity=self.verbosity))
        else:

            tree = TrClTree.cast(runPhyml(concat, self.tmpdir,
                                analysis=self.analysis,
                                verbosity=self.verbosity))
            if self.verbosity == 1:
                print

        # concat local variable dies here and goes to garbage collect

        self.concats[index_tuple] = tree
        return tree

    def concatenate(self, index_tuple):
        """ NB: had a version of this which used a reduce construct to
        concatenate the alignments - reduce(lambda x,y: x+y, records_list) - but
        this led to problems of the original object being modified. Deepcopying
        the first record, ensuring a new memory address for the concatenation,
        seems more robust. """

        member_records = self.members(index_tuple)
        concat = concatenate(member_records)
        concat.name = '-'.join(str(x) for x in index_tuple)
        return concat

    def update_history(self, score, index_tuple):
        time = timeit.default_timer()
        self.history.append([time, score, index_tuple, len(index_tuple)])

    def print_history(self, fh=sys.stdout):
        for iteration, (time, score, index_tuple, nclusters) in enumerate(
                self.history):
            fh.write(str(iteration) + "\t")
            fh.write(str(time) + "\t")
            fh.write(str(score) + "\t")
            fh.write(str(index_tuple) + "\t")
            fh.write(str(nclusters) + "\n")

    def clear_history(self):
        self.history = []

    def members(self, index_tuple):
        return [self.records[n] for n in index_tuple]

    def populate_cache(self):
        for i, rec in enumerate(self.records):
            key = (i,)
            tree = rec.tree
            self.concats[key]=tree

    def score(self, partition_object, history=True, **kwargs):
        """ Generates the index lists of the Partition object, gets the score
        for each one, and returns the sum """

        inds = partition_object.get_membership()
        self.add_partition_list([partition_object])
        likelihood = sum([self.add(index_tuple, **kwargs).score
                          for index_tuple in inds])
        if history is True:
            self.update_history(likelihood, inds)
        return(likelihood)

    def simulate(self, index_tuple, model=None):
        """ Simulate a group of sequence alignments using ALF. Uses one of
        {(GCB, JTT, LG, WAG - protein), (CPAM, ECM and ECMu - DNA)}, WAG by
        default. TO DO: add parameterised models when I have a robust (probably
        PAML) method of estimating them from alignment+tree """

        if self.datatype == 'protein':  # set some defaults
            model = model or 'WAG'
            optioncheck(model, [
                'CPAM',
                'ECM',
                'ECMu',
                'WAG',
                'JTT',
                'GCB',
                'LG',
                ])
        else:
            model = model or 'GTR'
            try:
                optioncheck(model, ['CPAM', 'ECM', 'ECMu', 'GTR'])
            except OptionError, e:
                print 'Choose a DNA-friendly model for simulation:\n', e
                return

        member_records = self.members(index_tuple)
        concat = self.concatenate(index_tuple)
        (lengths, names) = zip(*[(rec.seqlength, rec.name) for rec in
                               member_records])
        full_length = sum(lengths)
        concat.tree = self.add(index_tuple)

        simulated_records = simulate_from_record(
            concat,
            length=full_length,
            tmpdir=self.tmpdir,
            model=model,
            split_lengths=lengths,
            gene_names=names,
            )

        return simulated_records

    def simulate_from_result(self, partition_object, **kwargs):
        inds = partition_object.get_membership()
        return flatten_list([self.simulate(ind, **kwargs) for ind in inds])
