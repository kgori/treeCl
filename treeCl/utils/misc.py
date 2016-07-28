from __future__ import division
from builtins import str
from progressbar import ProgressBar, Percentage, SimpleProgress, Timer, AdaptiveETA, Bar, FormatLabel
import numpy as np
import itertools
import random
from phylo_utils import seq_to_partials
from phylo_utils.markov import TransitionMatrix
from phylo_utils.models import LG, WAG, GTR
from phylo_utils.likelihood import GammaMixture
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict

__all__ = ['concatenate',
           'flatten_list',
           'symmetrise',
           'regex_search_extract',
           'setup_progressbar',
           'model_translate',
           'smooth_freqs',
           'grouper',
           'insort_no_dup',
           'alignment_to_partials',
           'biopython_to_partials',
           'create_gamma_model',
           'weighted_choice',
           'sample_wr']

def concatenate(alignments):
    """
    Concatenates a list of Bio.Align.MultipleSeqAlignment objects.
    If any sequences are missing the are padded with unknown data
    (Bio.Seq.UnknownSeq).
    Returns a single Bio.Align.MultipleSeqAlignment.
    Limitations: any annotations in the sub-alignments are lost in
    the concatenated alignment.
    """

    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in alignments for seq in aln)

    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    tmp = defaultdict(list)

    # Assume all alignments have same alphabet
    alphabet = alignments[0]._alphabet

    for aln in alignments:
        length = aln.get_alignment_length()

        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels

        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the tmp dict
        for label in missing:
            new_seq = UnknownSeq(length, alphabet=alphabet)
            tmp[label].append(str(new_seq))

        # else stuff the string representation into the tmp dict
        for rec in aln:
            tmp[rec.id].append(str(rec.seq))

    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment
    msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(v), alphabet=alphabet), id=k, name=k, description=k)
               for (k,v) in tmp.items())
    return msa

def flatten_list(list_):
    newlist = list()
    x = newlist.extend
    ap = newlist.append
    for sublist in list_:
        try:
            x(sublist)
        except TypeError: # if the "sublist" is non-iterable, append as a plain element
            ap(sublist)
    return newlist


def symmetrise(matrix, tri='upper'):
    """
    Will copy the selected (upper or lower) triangle of a square matrix
    to the opposite side, so that the matrix is symmetrical.
    Alters in place.
    """
    if tri == 'upper':
        tri_fn = np.triu_indices
    else:
        tri_fn = np.tril_indices
    size = matrix.shape[0]
    matrix[tri_fn(size)[::-1]] = matrix[tri_fn(size)]
    return matrix


def regex_search_extract(search_attempt):
    return search_attempt.group() if search_attempt else None


def setup_progressbar(msg, size, format_label=None, simple_progress=False):
    if not msg.endswith(': '):
        msg += ': '

    if simple_progress:
        widgets = [msg,
                   SimpleProgress(), ' ',
                   Bar(), ' ',
                   Timer(), ' ',
                   AdaptiveETA()]
    else:
        widgets = [msg,
                   Percentage(), ' ',
                   Bar(), ' ',
                   Timer(), ' ',
                   AdaptiveETA()]

    if format_label is not None:
        widgets.append(FormatLabel(format_label))

    pbar = ProgressBar(widgets=widgets, maxval=size)
    return pbar

def model_translate(model):

    translation = {'LG' : 'LG08',
                   'WAG': 'WAG01'}
    return translation.get(model, model)

def smooth_freqs(freqs):
    """
    Smooths freqs vector, guarantees sum == 1
    :param freqs: vector of frequencies
    :return: vector of frequencies guaranteed to sum to 1
    """
    s = sum(freqs)
    return [f/s for f in freqs]

def grouper(n, iterable):
    """
    >>> list(grouper(3, 'ABCDEFG'))
    [['A', 'B', 'C'], ['D', 'E', 'F'], ['G']]
    """
    iterable = iter(iterable)
    return iter(lambda: list(itertools.islice(iterable, n)), [])

def insort_no_dup(lst, item):
    """
    If item is not in lst, add item to list at its sorted position
    """
    import bisect
    ix = bisect.bisect_left(lst, item)
    if lst[ix] != item: 
        lst[ix:ix] = [item]

def alignment_to_partials(alignment, missing_data=None):
    """ Generate a partials dictionary from a treeCl.Alignment """
    partials_dict = {}
    for (name, sequence) in alignment.get_sequences():
        datatype = 'dna' if alignment.is_dna() else 'protein'
        partials_dict[name] = seq_to_partials(sequence, datatype)

    if missing_data is not None:
        l = len(alignment)
        for name in missing_data:
            if name not in partials_dict:
                partials_dict[name] = seq_to_partials('-'*l, datatype)
    return partials_dict

def biopython_to_partials(alignment, datatype):
    """ Generate a partials dictionary from a treeCl.Alignment """
    partials_dict = {}
    for seq in alignment:
        partials_dict[seq.name] = seq_to_partials(seq, datatype)
    return partials_dict

def create_gamma_model(alignment, missing_data=None, ncat=4):
    """ Create a phylo_utils.likelihood.GammaMixture for calculating
    likelihood on a tree, from a treeCl.Alignment and its matching 
    treeCl.Parameters """
    model = alignment.parameters.partitions.model
    freqs = alignment.parameters.partitions.frequencies
    alpha = alignment.parameters.partitions.alpha
    if model == 'LG':
        subs_model = LG(freqs)
    elif model == 'WAG':
        subs_model = WAG(freqs)
    elif model == 'GTR':
        rates = alignment.parameters.partitions.rates
        subs_model = GTR(rates, freqs, True)
    else:
        raise ValueError("Can't handle this model: {}".format(model))
    tm = TransitionMatrix(subs_model)
    gamma = GammaMixture(alpha, ncat)
    gamma.init_models(tm, alignment_to_partials(alignment, missing_data))
    return gamma

def weighted_choice(choices):
    total = sum(w for c, w in choices)
    r = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        if upto + w > r:
            return c
        upto += w
    assert False, "Shouldn't get here"

def sample_wr(lst):
    """
    Sample from lst, with replacement
    """
    arr = np.array(lst)
    indices = np.random.randint(len(lst), size=len(lst))
    sample = np.empty(arr.shape, dtype=arr.dtype)
    for i, ix in enumerate(indices):
        sample[i] = arr[ix]
    return list(sample)
