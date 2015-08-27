import collections
import os
import tempfile

from numpy import log

import bpp
from treeCl.utils import pll_helpers
from parameters import Parameters, PartitionParameters
from utils import fileIO
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from collections import defaultdict


class Alignment(bpp.Alignment):
    def __init__(self, *args):
        super(Alignment, self).__init__(*args)
        self.infile = None
        self.name = None
        self.parameters = Parameters()
        if len(args) > 0 and isinstance(args[0], basestring) and fileIO.can_locate(args[0]):
            self.infile = os.path.abspath(args[0])
            self.parameters.filename = args[0]
            self.name = os.path.splitext(os.path.basename(self.infile))[0]

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

    def write_alignment(self, filename, file_format, interleaved=None):
        """
        Overloads the base class function.
        Uses Bio AlignIO to write because biopp writes
        phylip interleaved in a way that causes an error
        with FastTree
        """
        b = self.to_biopython_msa()
        if file_format == 'phylip':
            file_format = 'phylip-relaxed'
        AlignIO.write(b, filename, file_format)

    def get_alignment_file(self, as_phylip=True):
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
            return pll_helpers.create_instance(alignment, *args)  # args=(partitions, tree, threads, rns)

        except (IOError, TypeError):
            with fileIO.TempFile() as tmpfile:
                self.write_alignment(tmpfile, 'phylip', True)
                return pll_helpers.create_instance(tmpfile, *args)

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

    def to_biopython_msa(self):
        alph = IUPAC.extended_dna if self.is_dna() else IUPAC.extended_protein
        msa = MultipleSeqAlignment([SeqRecord(Seq(sequence, alphabet=IUPAC.extended_dna), id=key) for (key, sequence) in self.get_sequences()])
        for seq in msa: seq.description=''
        return msa

