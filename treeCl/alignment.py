import collections
import itertools
import os
import tempfile

import numpy as np

import bpp
from parameters import Parameters, PartitionParameters
from utils import fileIO, alignment_to_partials
from distance_matrix import DistanceMatrix
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from collections import defaultdict
from phylo_utils.likcalc import discrete_gamma, sitewise_lik, sitewise_lik_derivs
from phylo_utils import seq_to_partials
from phylo_utils.markov import TransitionMatrix
from phylo_utils.models import LG, GTR
import logging
logger = logging.getLogger(__name__)

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

    def get_unconstrained_likelihood(self):
        weights = collections.Counter()
        sites = self.get_sites()
        n = len(sites)
        ucl = 0

        for site in sites:
            weights[site] += 1

        for v in weights.itervalues():
            if v > 1:
                ucl += v * np.log(v)
        return ucl - n * np.log(n)

    def to_biopython_msa(self):
        alph = IUPAC.extended_dna if self.is_dna() else IUPAC.extended_protein
        msa = MultipleSeqAlignment([SeqRecord(Seq(sequence, alphabet=IUPAC.extended_dna), id=key) for (key, sequence) in self.get_sequences()])
        for seq in msa: seq.description=''
        return msa

def pairdists(alignment, ncat=4, tolerance=1e-6):
    """ Load an alignment, calculate all pairwise distances and variances """
    def calc(brlen):
        """
        Inner function calculates l'hood + derivs at branch length = brlen
        """
        result = sum([sitewise_lik_derivs(tm.get_p_matrix(gamma_rates[k]*brlen),
                                          tm.get_dp_matrix(gamma_rates[k]*brlen),
                                          tm.get_d2p_matrix(gamma_rates[k]*brlen),
                                          tm.freqs, partials[key1], partials[key2])*(1.0/ncat)
                              for k in range(ncat)])
        lk = np.log(result[:,0]).sum()
        dlk = (result[:,1]/result[:,0]).sum()
        d2lk = ((result[:,0]*result[:,2] - result[:,1]**2)/result[:,0]**2).sum()
        return lk, dlk, d2lk

    def get_step(dlk, d2lk):
        step = dlk / np.abs(d2lk) # abs makes optimiser backtrack from a minimum likelihood
        while (brlen + step) < 0:
            step *= 0.5
        return step

    try:
        model = alignment.parameters.partitions.model
        freqs = alignment.parameters.partitions.frequencies
        alpha = alignment.parameters.partitions.alpha
    except:
        logger.error('No parameters available')
        return

    if model == 'LG':
        subs_model = LG(freqs)
    elif model == 'GTR':
        rates = alignment.parameters.partitions.rates
        subs_model = GTR(rates, freqs, True)
    else:
        raise ValueError("Can't handle this model: {}".format(model))

    # Set up markov model
    tm = TransitionMatrix(subs_model)
    gamma_rates = discrete_gamma(alpha, ncat)
    partials = alignment_to_partials(alignment)
    seqnames = alignment.get_names()
    nseq = len(seqnames)
    distances = np.zeros((nseq, nseq))
    variances = np.zeros((nseq, nseq))

    for i, j in itertools.combinations(range(nseq), 2):
        key1 = seqnames[i]
        key2 = seqnames[j]
        maxiter = 100

        brlen = 1.0  # just have a guess
        lk, dlk, d2lk = calc(brlen)
        maxlk = lk
        niter = 0
        step = get_step(dlk, d2lk)

        # This is the newton optimiser
        while True:
            niter += 1
            if niter > maxiter:
                break  # failed to converge somehow

            # Do the calculation to work out the new step
            lk, dlk, d2lk = calc(brlen + step)
            if (lk - maxlk) < -1000*tolerance:
                # the likelihood got worse, so the step was too big
                # so restore the old values and halve the step, try again
                step *= 0.5
                continue

            else:
                # successful move. update brlen
                brlen = brlen + step
                maxlk = lk
                step = get_step(dlk, d2lk)

            if np.abs(dlk) < tolerance:
                break  # Converged

        distances[i, j] = distances[j, i] = brlen
        variances[i, j] = variances[j, i] = np.abs(-1/d2lk)
    dm = DistanceMatrix.from_array(distances, names=seqnames)
    vm = DistanceMatrix.from_array(variances, names=seqnames)
    return dm, vm
