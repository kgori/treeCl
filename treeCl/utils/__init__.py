# import fileIO
from progressbar import ProgressBar, Percentage, Timer, AdaptiveETA, Bar, FormatLabel
import numpy as np

from printing import print_and_return


def flatten_list(list_):
    newlist = list()
    x = newlist.extend
    for sublist in list_:
        x(sublist)
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
    return (search_attempt.group() if search_attempt else None)


def setup_progressbar(msg, size, format_label=None):
    if not msg.endswith(': '):
        msg += ': '

    widgets = [msg,
               Percentage(), ' ',
               Bar(), ' ',
               Timer(), ' ',
               AdaptiveETA()]

    if format_label is not None:
        widgets.append(FormatLabel(format_label))

    pbar = ProgressBar(widgets=widgets, maxval=size)
    return pbar
