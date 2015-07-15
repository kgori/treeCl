# import fileIO
from progressbar import ProgressBar, Percentage, SimpleProgress, Timer, AdaptiveETA, Bar, FormatLabel
import numpy as np

from printing import print_and_return
from Bio.Seq import Seq, UnknownSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict

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
    msa = MultipleSeqAlignment(SeqRecord(Seq(''.join(v), alphabet=alphabet), id=k)
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
    return translation[model]

def smooth_freqs(freqs):
    """
    Smooths freqs vector, guarantees sum == 1
    :param freqs: vector of frequencies
    :return: vector of frequencies guaranteed to sum to 1
    """
    s = sum(freqs)
    return [f/s for f in freqs]
