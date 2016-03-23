from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from builtins import object
import collections
import itertools
import os

import numpy as np
from six import string_types

from .parameters import Parameters
from .utils import fileIO, alignment_to_partials, concatenate, sample_wr
from .distance_matrix import DistanceMatrix
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from phylo_utils.likcalc import discrete_gamma, sitewise_lik, sitewise_lik_derivs
from phylo_utils.markov import TransitionMatrix
from phylo_utils.models import LG, GTR
from phylo_utils.likcalc import _evolve_states, _weighted_choices
import logging
logger = logging.getLogger(__name__)

def set_alphabet(msa, alphabet):
    msa._alphabet = alphabet
    for seqrec in msa:
        seqrec.seq.alphabet = alphabet


class Alignment(object):
    def __init__(self, *args, **kwargs):
        """
        Initialise an alignment.
        Alignment()                                - initialise an empty alignment
        Alignment([alignments...])                 - concatenate a list of alignments into one alignment
        Alignment([('seq1', 'acgt...'),            - construct from a list of (header, sequence) tuples
                   ('seq2', 'tgca...'),
                   ...])
        Alignment(file_path, (file_format))        - read an alignment from file
        Alignment(..., alphabet=<'dna'|'protein'>) - specify the alphabet, either dna or protein
        """
        self.infile = None
        self.name = None
        self.parameters = Parameters()
        if len(args) == 0:
            self._msa = None
        elif isinstance(args[0], list):
            # initialise from a list of sequences or Alignments
            if all(isinstance(element, self.__class__) for element in args[0]):
                # initialise from list of Alignment objects
                self._msa = concatenate([al._msa for al in args[0]])
            elif all(isinstance(element, tuple) for element in args[0]):
                # initialise from list of (name, sequence) tuples
                msa = MultipleSeqAlignment([SeqRecord(Seq(sequence), id=key, description=key, name=key) 
                                            for (key, sequence) in args[0]])
                self._msa = self._guess_alphabet(msa)
                if 'name' in kwargs:
                    self.name = kwargs['name']
            else:
                self._msa = None

        elif len(args) > 0 and isinstance(args[0], string_types) and fileIO.can_locate(args[0]):
            self.infile = os.path.abspath(args[0])
            self.parameters.filename = args[0]
            self.name = os.path.splitext(os.path.basename(self.infile))[0]
            if args[1]=='phylip':
                if len(args) >= 3:
                    if args[2]:
                        self.read_alignment(args[0], 'phylip-relaxed')
                    else:
                        self.read_alignment(args[0], 'phylip-sequential')
                self.read_alignment(args[0], 'phylip-relaxed')
            else:
                self.read_alignment(args[0], args[1])

        else:
            logger.warn('Failed to initialise alignment - couldn\'t read args as a file or interpret as sequences')
            self._msa = None

        if 'alphabet' in kwargs and self._msa is not None:
            alphabet = kwargs['alphabet']
            if alphabet in ('dna', 'DNA'):
                set_alphabet(self._msa, IUPAC.ambiguous_dna)
            elif alphabet in ('protein', 'PROTEIN'):
                set_alphabet(self._msa, IUPAC.extended_protein)
            else:
                logger.warning('Set alphabet to "dna" or "protein", not {}'.format(alphabet))

    def __add__(self, other):
        return self.__class__([self, other])

    def __str__(self):
        contents = self.get_sequences()
        output_string = 'Alignment: {}\n'.format(self.name)

        return output_string + ''.join(
            ['>{}\n{}\n... ({}) ...\n{}\n'.format(header, seq[:50], len(seq) - 100, seq[-50:]) for header, seq in
             contents])

    def __len__(self):
        if self._msa:
            return self._msa.get_alignment_length()

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

    def is_dna(self):
        return isinstance(self._msa._alphabet, (type(IUPAC.ambiguous_dna), type(IUPAC.unambiguous_dna)))

    def is_protein(self):
        return isinstance(self._msa._alphabet, (type(IUPAC.protein), type(IUPAC.extended_protein)))
    
    def read_alignment(self, *args, **kwargs):
        filename = args[0]
        args = args[1:]
        with open(filename) as fl:
            msa = AlignIO.read(fl, *args, **kwargs)
        self.infile = filename
        # guess alphabet
        self._msa = self._guess_alphabet(msa)

    def _guess_alphabet(self, msa):
        allchars = [char.upper() for sr in msa for char in str(sr.seq) if char.upper() not in '-?X']
        probably_dna = set(allchars).issubset(set(IUPAC.ambiguous_dna.letters))
        if probably_dna:
            alphabet = IUPAC.ambiguous_dna
        else:
            alphabet = IUPAC.extended_protein
        set_alphabet(msa, alphabet)
        return msa

    def write_alignment(self, filename, file_format, interleaved=None):
        """
        Write the alignment to file using Bio.AlignIO
        """
        if file_format == 'phylip':
            file_format = 'phylip-relaxed'
        AlignIO.write(self._msa, filename, file_format)

    def get_alignment_file(self, as_phylip=True):
        try:
            with open(self.infile) as fl:
                if as_phylip:
                    header = fl.readline().strip().split()
                    assert len(header) == 2 and header[0].isdigit() and header[1].isdigit()
            alignment = self.infile
            return os.path.abspath(alignment), False

        except (IOError, TypeError, AssertionError):
            with fileIO.NonDeletingTempFile() as tmpfile:
                self.write_alignment(tmpfile, "phylip", interleaved=True)
            return tmpfile, True

    def get_unconstrained_likelihood(self):
        weights = collections.Counter()
        sites = self.get_sites()
        n = len(sites)
        ucl = 0

        for site in sites:
            weights[site] += 1

        for v in weights.values():
            if v > 1:
                ucl += v * np.log(v)
        return ucl - n * np.log(n)

    def to_biopython_msa(self):
        return self._msa
        
    def get_sequences(self):
        if self._msa:
            return [(sr.name, str(sr.seq)) for sr in self._msa]

    def get_sites(self):
        if self._msa:
            seqs = [str(sr.seq) for sr in self._msa]
            return [''.join(col) for col in zip(*seqs)]

    def get_names(self):
        if self._msa:
            return [sr.name for sr in self._msa]

    def compute_distances(self, ncat=4, tolerance=1e-6):
        try:
            return pairdists(self, ncat, tolerance)
        except:
            return

    def simulate(self, nsites, transition_matrix, tree, ncat=1, alpha=1):
        """
        Return sequences simulated under the transition matrix's model 
        """
        sim = SequenceSimulator(transition_matrix, tree, ncat, alpha)
        return list(sim.simulate(nsites).items())

    def bootstrap(self):
        """
        Return a new Alignment that is a bootstrap replicate of self
        """
        new_sites = sorted(sample_wr(self.get_sites()))
        seqs = list(zip(self.get_names(), (''.join(seq) for seq in zip(*new_sites))))
        return self.__class__(seqs)


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
        variances[i, j] = variances[j, i] = np.abs(-1.0/d2lk)
    dm = DistanceMatrix.from_array(distances, names=seqnames)
    vm = DistanceMatrix.from_array(variances, names=seqnames)
    return dm, vm


class SequenceSimulator(object):
    """
    Simple markov generator to produce tip sequences on a tree
    according to a transition matrix
    """
    def __init__(self, transmat, tree, ncat=1, alpha=1):
        """
        Initialise the simulator with a transition matrix and a tree.
        The tree should have branch lengths. If it doesn't this will
        trigger a warning, but will continue.
        """
        # store the tree
        self.tree = tree
        self.states = np.array(transmat.model.states)
        self.state_indices = np.array(list(range(transmat.model.size)), dtype=np.intc)
        # initialise equilibrium frequency distribution
        self.freqs = transmat.freqs
        # Gamma rate categories
        self.ncat = ncat
        self.alpha = alpha
        self.gamma_rates = discrete_gamma(alpha, ncat)
        
        # initialise probabilities on tree
        for node in self.tree.preorder(skip_seed=True):
            l = node.edge.length or 0
            if l == 0:
                print ('warning')
                #logger.warn('This tree has zero length edges')
            nstates = self.states.shape[0]
            node.pmats = np.empty((ncat, nstates, nstates))
            for i in range(ncat):
                node.pmats[i] = transmat.get_p_matrix(l*self.gamma_rates[i])

        self.sequences = {}

    def simulate(self, n):
        """
        Evolve multiple sites during one tree traversal
        """
        self.tree._tree.seed_node.states = self.ancestral_states(n)
        categories = np.random.randint(self.ncat, size=n).astype(np.intc)

        for node in self.tree.preorder(skip_seed=True):
            node.states = self.evolve_states(node.parent_node.states, categories, node.pmats)
            if node.is_leaf():
                self.sequences[node.taxon.label] = node.states
        return self.sequences_to_string()

    def ancestral_states(self, n):
        """
        Generate ancestral sequence states from the equilibrium frequencies
        """
        anc = np.empty(n, dtype=np.intc)
        _weighted_choices(self.state_indices, self.freqs, anc)
        return anc

    def evolve_states(self, parent_states, categories, probs):
        """
        Evolve states from parent to child. States are sampled
        from gamma categories passed in the array 'categories'.
        The branch length information is encoded in the probability
        matrix, 'probs', generated in __init__.
        """
        child_states = np.empty(parent_states.shape, dtype=np.intc)
        _evolve_states(self.state_indices, parent_states, categories, probs, child_states)
        return child_states
    
    def sequences_to_string(self):
        """
        Convert state indices to a string of characters
        """
        return {k: ''.join(self.states[v]) for (k, v) in self.sequences.items()}
