import bpp
import collections
import os
import shutil
import sys
import tempfile

from numpy import log

from interfacing import pll
from parameters import Parameters
from utils import fileIO


class Alignment(bpp.Alignment):
    def __init__(self, *args):
        super(Alignment, self).__init__(*args)
        self.infile = None
        self.name = None
        self.parameters = Parameters()
        if len(args) > 0 and isinstance(args[0], basestring) and fileIO.can_locate(args[0]):
            self.infile = args[0]
            self.parameters.filename = args[0]

    def __add__(self, other):
        return self.__class__([self, other])

    def __str__(self):
        contents = self.get_sequences()
        output_string = 'Alignment: {}\n'.format(self.name)

        return output_string + ''.join(
            ['>{}\n{}\n... ({}) ...\n{}\n'.format(header, seq[:50], len(seq) - 100, seq[-50:]) for header, seq in
             contents])

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, parameters):
        self._parameters = parameters

    @property
    def tree(self):
        if self.parameters.ml_tree is not None:
            return self.parameters.ml_tree
        elif self.parameters.nj_tree is not None:
            return self.parameters.nj_tree
        else:
            raise AttributeError('No tree')

    def read_alignment(self, *args, **kwargs):
        super(Alignment, self).read_alignment(*args, **kwargs)
        self.infile = args[0]

    def get_alignment_file(self, as_phylip=False):
        try:
            with open(self.infile) as fl:
                if as_phylip:
                    header = fl.readline().strip().split()
                    assert len(header) == 2 and header[0].isdigit() and header[1].isdigit()
            alignment = self.infile
            return os.path.abspath(alignment), False

        except (IOError, TypeError, AssertionError):
            tmpfile = os.path.abspath(tempfile.mkstemp()[1])
            self.write_alignment(tmpfile, "phylip", interleaved=True)
            return tmpfile, True

    def pll_get_instance(self, *args):
        try:
            with open(self.infile):
                pass
            alignment = self.infile
            return pll.create_instance(alignment, *args)  # args=(partitions, tree, threads, rns)

        except (IOError, TypeError):
            with fileIO.TempFile() as tmpfile:
                self.write_alignment(tmpfile, 'phylip', True)
                return pll.create_instance(tmpfile, *args)

    def get_unconstrained_likelihood(self):
        weights = collections.Counter()
        sites = self.get_sites()
        n = len(sites)
        ucl = 0

        for site in sites:
            weights[site] += 1

        for v in weights.itervalues():
            if v > 1:
                ucl += v * log(v)
        return ucl - n * log(n)

    def set_params_from_pll_result(self, result):
        self.parameters.ml_tree = result['ml_tree']
        self.parameters.likelihood = result['likelihood']
        params = PartitionParameters()
        params.alpha = result['partitions'][0]['alpha']
        params.frequencies = result['partitions'][0]['frequencies']
        try:
            params.rates = result['partitions'][0]['rates']
        except KeyError:
            pass
        params.name = result['partitions'][0]['name']
        self.parameters.partitions = [params]
