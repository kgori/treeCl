#!/usr/bin/env python

import glob
import re
from copy import deepcopy
from treeCl.sequence_record import TCSeqRec, concatenate
from treeCl.externals import runDV, runTC, runPhyml, simulate_from_tree
from treeCl.distance_matrix import DistanceMatrix
from errors import FileError, DirectoryError, OptionError, \
    optioncheck_and_raise, directorycheck_and_make, directorycheck_and_raise
from utils import fileIO

sort_key = lambda item: tuple((int(num) if num else alpha) for (num, alpha) in
                              re.findall(r'(\d+)|(\D+)', item))

get_name = lambda i: i[i.rindex('/') + 1:i.rindex('.')]


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
        tmpdir='/tmp',
        calc_distances=False,
        ):

        self.tmpdir = directorycheck_and_make(tmpdir)

        if records:
            self.records = records
            self.datatype = datatype or records[0].datatype
            optioncheck_and_raise(self.datatype, ['dna', 'protein'])

        elif input_dir:
            directorycheck_and_raise(input_dir)
            self.datatype = optioncheck_and_raise(datatype, ['dna', 'protein'])
            optioncheck_and_raise(file_format, ['fasta', 'phylip'])
            self.records = self.read_files(input_dir, file_format)

        else:
            print 'Provide a list of records, or the path to a set of alignments'

        if calc_distances:
            self.calc_distances()

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

    def read_files(self, input_dir, file_format):
        """ Get list of alignment files from an input directory *.fa, *.fas and
        *.phy files only
        
        Stores in self.files """

        if file_format == 'fasta':
            files = glob.glob('{0}/*.fa'.format(input_dir))
            if len(files) == 0:
                files = glob.glob('{0}/*.fas'.format(input_dir))
        elif file_format == 'phylip':
            files = glob.glob('{0}/*.phy'.format(input_dir))
        else:
            print 'Unrecognised file format %s' % file_format
            files = None
        if not files:
            print 'No sequence files found in {0}'.format(input_dir)
            raise FileError(input_dir)
        files = sorted(files, key=sort_key)
        return [TCSeqRec(f, file_format=file_format, datatype=self.datatype,
                name=get_name(f), tmpdir=self.tmpdir) for f in files]

    def calc_distances(self, verbosity=0):
        for rec in self.records:
            runDV(rec, verbosity)

    def calc_TC_trees(self, verbosity=0):
        for rec in self.records:
            runTC(rec, verbosity=verbosity)

    def calc_ML_trees(self, verbosity=0):
        for rec in self.records:
            runPhyml(rec, analysis='ml', verbosity=verbosity)

    def calc_NJ_trees(self, verbosity=0):
        for rec in self.records:
            runPhyml(rec, analysis='nj', verbosity=verbosity)

    def distance_matrix(self, metric):
        """ Generate a distance matrix from a fully-populated Collection """

        trees = [rec.tree for rec in self.records]
        return DistanceMatrix(trees, metric, tmpdir=self.tmpdir)

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
        max_guidetrees=10,
        tmpdir=None,
        datatype=None,
        ):

        self.analysis = optioncheck_and_raise(analysis, ['ml', 'nj',
                'TreeCollection'])
        self.max_guidetrees = max_guidetrees
        self.records = records
        self.datatype = datatype or records[0].datatype
        optioncheck_and_raise(self.datatype, ['protein', 'dna'])
        self.tmpdir = tmpdir or records[0].tmpdir
        self.concats = {}

    def add(self, index_list, verbosity=1):
        """ Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree """

        if index_list in self.concats:
            return self.concats[index_list]

        concat = self.concatenate(index_list)

        if self.analysis == 'TreeCollection':
            guidetrees = [self.records[n].tree for n in
                          index_list][:self.max_guidetrees]
            tree = runTC(concat, guidetrees, verbosity=verbosity)
        else:

            tree = runPhyml(concat, analysis=self.analysis, verbosity=verbosity)

        # concat local variable dies here and goes to garbage collect

        self.concats[index_list] = tree
        return tree

    def concatenate(self, index_list):
        """ NB: had a version of this which used a reduce construct to
        concatenate the alignments - reduce(lambda x,y: x+y, records_list) - but
        this led to problems of the original object being modified. Deepcopying
        the first record, ensuring a new memory address for the concatenation,
        seems more robust. """

        member_records = self.members(index_list)
        concat = concatenate(member_records)
        concat.name = '-'.join(str(x) for x in index_list)
        return concat

    def members(self, index_list):
        return [self.records[n] for n in index_list]

    def score(self, partition_object, **kwargs):
        """ Generates the index lists of the Partition object, gets the score
        for each one, and returns the sum """

        inds = partition_object.get_membership()
        return sum([self.add(ind, **kwargs).score for ind in inds])

    def simulate(self, index_list, model=None):
        """ Simulate a group of sequence alignments using ALF. Uses one of
        {(GCB, JTT, LG, WAG - protein), (CPAM, ECM and ECMu - DNA)}, WAG by
        default. TO DO: add parameterised models when I have a robust (probably
        PAML) method of estimating them from alignment+tree """

        if self.datatype == 'protein':  # set some defaults
            model = model or 'WAG'
            optioncheck_and_raise(model, [
                'CPAM',
                'ECM',
                'ECMu',
                'WAG',
                'JTT',
                'GCB',
                'LG',
                ])
        else:
            model = model or 'ECM'
            try:
                optioncheck_and_raise(model, ['CPAM', 'ECM', 'ECMu'])
            except OptionError, e:
                print 'Choose a DNA-friendly model for simulation:\n', e
                return

        member_records = self.members(index_list)
        (lengths, names) = zip(*[(rec.seqlength, rec.name) for rec in
                               member_records])
        full_length = sum(lengths)
        tree = self.add(index_list)

        simulated_records = simulate_from_tree(
            tree=tree,
            length=full_length,
            datatype=self.datatype,
            tmpdir=self.tmpdir,
            model=model,
            split_lengths=lengths,
            gene_names=names,
            )

        return simulated_records

    def simulate_from_result(self, partition_object, **kwargs):
        inds = partition_object.get_membership()
        return [self.simulate(ind, **kwargs) for ind in inds]
