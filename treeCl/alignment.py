from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import str
from builtins import zip
from builtins import range
from builtins import object
import collections
import itertools
import io, os, random

import numpy as np
import pandas as pd
from six import string_types

from .parameters import Parameters
from .constants import ISPY3
from .utils import fileIO, alignment_to_partials, concatenate, sample_wr
from .distance_matrix import DistanceMatrix
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC
from Bio import AlignIO
from phylo_utils.likcalc import discrete_gamma
from phylo_utils.markov import TransitionMatrix
import phylo_utils.models
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
            logger.warning('Failed to initialise alignment - couldn\'t read args as a file or interpret as sequences')
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
        with fileIO.freader(filename) as fl:
            if ISPY3:
                handle = io.TextIOWrapper(fl)
                msa = AlignIO.read(handle, *args, **kwargs)
            else:
                msa = AlignIO.read(fl, *args, **kwargs)
        self.infile = filename
        # guess alphabet
        self._msa = self._guess_alphabet(msa)

    def _guess_alphabet(self, msa):
        if msa.get_alignment_length() > 1000:
            allchars = [char for sr in msa for char in random.sample(list(sr.seq.upper()), 1000)]
        elif len(msa) > 100:
            allchars = [char for sr in random.sample(list(msa), 100) for char in sr.seq.upper()]
        elif len(msa) > 100 and msa.get_alignment_length() > 1000:
            allchars = [char for sr in random.sample(list(msa), 100)
                            for char in random.sample(list(sr.seq.upper(), 1000))]
        else:
            allchars = [char for sr in msa for char in str(sr.seq.upper())]
        probably_dna = (set(allchars) - set('-?X')).issubset(set(IUPAC.ambiguous_dna.letters))
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
            with fileIO.TempFile(disable_delete=True) as tmpfile:
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

    def compute_distances(self, model, alpha=None, ncat=4, tolerance=1e-6):
        """
        Compute pairwise distances between all sequences according to
        a given substitution model, `model`, of type phylo_utils.models.Model
        (e.g. phylo_utils.models.WAG(freqs), phylo_utils.models.GTR(rates, freqs))
        The number of gamma categories is controlled by `ncat`. Setting ncat=1
        disable gamma rate variation. The gamma alpha parameter must be supplied
        to enable gamma rate variation.
        """
        return pairdists(self, model, alpha, ncat, tolerance)


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

    def to_data_frame(self):
        return pd.DataFrame([list(sr.seq) for sr in self._msa], index=self.get_names())

    @classmethod
    def from_data_frame(cls, df):
        return cls([(label, ''.join(row)) for (label, row) in df.iterrows()])


class BranchLengthOptimiser(object):
    """
    Wrapper for use with scipy optimiser (e.g. brenth/brentq)
    """

    def __init__(self, node1, node2, initial_brlen=1.0):
        self.root = node1
        self.desc = node2
        self.updated = None
        self.__call__(initial_brlen)

    def __call__(self, brlen):
        if brlen < 0:
            self.lnl, self.dlnl, self.d2lnl = -np.inf, np.nan, np.nan
            self.updated = brlen
            return self.lnl, self.dlnl, self.d2lnl
        if self.updated != brlen:
            self.updated = brlen
            self.lnl, self.dlnl, self.d2lnl = self.root.compute_likelihood(self.desc, brlen, derivatives=True)
        return self.lnl, self.dlnl, self.d2lnl

    def get_lnl(self, brlen):
        return self.__call__(brlen)[0]

    def get_dlnl(self, brlen):
        return np.array([self.__call__(brlen)[1]])

    def get_d2lnl(self, brlen):
        return np.array([self.__call__(brlen)[2]])

    def get_negative_lnl(self, brlen):
        return -self.__call__(max(0,brlen))[0]

    def get_negative_dlnl(self, brlen):
        return -self.__call__(max(0,brlen))[1]

    def get_negative_d2lnl(self, brlen):
        return -self.__call__(max(0,brlen))[2]

    def __str__(self):
        return 'Branch length={}, Variance={}, Likelihood+derivatives = {} {} {}'.format(self.updated, -1 / self.d2lnl,
                                                                                         self.lnl, self.dlnl,
                                                                                         self.d2lnl)

def brent_optimise(node1, node2, min_brlen=0.001, max_brlen=10, verbose=False):
    """
    Optimise ML distance between two partials. min and max set brackets
    """
    from scipy.optimize import minimize_scalar
    wrapper = BranchLengthOptimiser(node1, node2, (min_brlen + max_brlen) / 2.)
    n = minimize_scalar(lambda x: -wrapper(x)[0], method='brent', bracket=(min_brlen, max_brlen))['x']
    if verbose:
        logger.info(wrapper)
    if n < min_brlen:
        n = min_brlen
        wrapper(n)
    return n, -1 / wrapper.get_d2lnl(n)

def pairdists(alignment, subs_model, alpha=None, ncat=4, tolerance=1e-6, verbose=False):
    """ Load an alignment, calculate all pairwise distances and variances
        model parameter must be a Substitution model type from phylo_utils """

    # Check
    if not isinstance(subs_model, phylo_utils.models.Model):
        raise ValueError("Can't handle this model: {}".format(model))

    if alpha is None:
        alpha = 1.0
        ncat = 1

    # Set up markov model
    tm = TransitionMatrix(subs_model)

    gamma_rates = discrete_gamma(alpha, ncat)
    partials = alignment_to_partials(alignment)
    seqnames = alignment.get_names()
    nseq = len(seqnames)
    distances = np.zeros((nseq, nseq))
    variances = np.zeros((nseq, nseq))

    # Check the model has the appropriate size
    if not subs_model.size == partials[seqnames[0]].shape[1]:
        raise ValueError("Model {} expects {} states, but the alignment has {}".format(model.name,
                                                                                       model.size,
                                                                                       partials[seqnames[0]].shape[1]))

    nodes = [phylo_utils.likelihood.LnlModel(tm) for seq in range(nseq)]
    for node, header in zip(nodes, seqnames):
        node.set_partials(partials[header])  # retrieve partial likelihoods from partials dictionary

    for i, j in itertools.combinations(range(nseq), 2):
        brlen, var = brent_optimise(nodes[i], nodes[j], verbose=verbose)
        distances[i, j] = distances[j, i] = brlen
        variances[i, j] = variances[j, i] = var
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
                logger.warning('This tree has zero length edges')
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
