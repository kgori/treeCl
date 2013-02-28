#!/usr/bin/env python

import glob
import re
from copy import deepcopy
from treeCl.sequence_record import TCSeqRec
from treeCl.externals import runDV, runTC, runPhyml
from treeCl.distance_matrix import DistanceMatrix
from errors import FileError, DirectoryError, optioncheck_and_raise
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
        input_dir,
        file_format,
        datatype,
        tmpdir='/tmp',
        calc_distances=False,
        ):

        if not fileIO.can_open(input_dir):  # checks
            raise DirectoryError(input_dir)
        if not fileIO.can_open(tmpdir):
            raise DirectoryError(tmpdir)

        self.file_format = file_format
        self.datatype = datatype
        self.tmpdir = tmpdir

        # files = self.read_files()

        self.records = input_dir
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
    def records(self, input_dir):
        files = self.read_files(input_dir)
        records = [TCSeqRec(f, file_format=self.file_format,
                   datatype=self.datatype, name=get_name(f),
                   tmpdir=self.tmpdir) for f in files]
        for rec in records:
            rec.sanitise()
        self._records = dict(enumerate(records))

        # self.reverse_lookup = {v:k for (k,v) in self.records}

    def read_files(self, input_dir):
        """ Get list of alignment files from an input directory *.fa, *.fas and
        *.phy files only
        
        Stores in self.files """

        if self.file_format == 'fasta':
            files = glob.glob('{0}/*.fa'.format(input_dir))
            if len(files) == 0:
                files = glob.glob('{0}/*.fas'.format(input_dir))
        elif self.file_format == 'phylip':
            files = glob.glob('{0}/*.phy'.format(input_dir))
        else:
            print 'Unrecognised file format %s' % self.file_format
            files = None
        if not files:
            print 'No sequence files found in {0}'.format(input_dir)
            raise FileError(input_dir)
        return sorted(files, key=sort_key)

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


class Scorer(object):

    """ Takes an index list, generates a concatenated SequenceRecord, calculates
    a tree and score """

    def __init__(
        self,
        records,
        analysis,
        max_guidetrees=10,
        ):
        optioncheck_and_raise(analysis, ['ml', 'nj', 'TreeCollection'])
        self.concats = {}
        self.records = records
        self.analysis = analysis
        self.max_guidetrees = max_guidetrees

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

        member_records = [self.records[n] for n in index_list]
        first_rec = member_records.pop(0)
        seed = deepcopy(first_rec)  # use of deepcopy here
        seed.tmpdir = first_rec.tmpdir
        for rec in member_records:  # is important
            seed += rec
        seed.name = '-'.join(str(x) for x in index_list)
        return seed

    def score(self, partition_object, **kwargs):
        """ Generates the index lists of the Partition object, gets the score
        for each one, and returns the sum """

        inds = partition_object.get_membership()
        return sum([self.add(ind, **kwargs).score for ind in inds])
