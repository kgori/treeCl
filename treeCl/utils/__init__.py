# import fileIO
from progressbar import ProgressBar, Percentage, SimpleProgress, Timer, AdaptiveETA, Bar, FormatLabel
import numpy as np

from printing import print_and_return


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
