#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div

# standard lib
import sys

# third party
import numpy as np
import pandas as pd
import sklearn

# treeCl
from . import errors
from .utils import fileIO

# Utilities
def isconnected(mask):
    """ Checks that all nodes are reachable from the first node - i.e. that the
    graph is fully connected. """

    nodes_to_check = list((np.where(mask[0, :])[0])[1:])
    seen = [True] + [False] * (len(mask) - 1)
    while nodes_to_check and not all(seen):
        node = nodes_to_check.pop()
        reachable = np.where(mask[node, :])[0]
        for i in reachable:
            if not seen[i]:
                nodes_to_check.append(i)
                seen[i] = True
    return all(seen)


def affinity(
        matrix,
        mask=None,
        scale=None,
):
    """
    Mask is a 2d boolean matrix. Scale is a 2d local scale matrix,
    as output by kscale(). It's the outer product of the kdists column
    vector produced by kdists.
    """

    mask = (mask if mask is not None else np.ones(matrix.shape, dtype=bool))
    assert isconnected(mask)
    scale = (scale if scale is not None else np.ones(matrix.shape))

    ix = np.where(np.logical_not(mask))
    scaled_matrix = -matrix ** 2 / scale
    # inputs where distance = 0 and scale = 0 result in NaN:
    # the next line replaces NaNs with -1.0
    scaled_matrix[np.where(np.isnan(scaled_matrix))] = -1.0
    affinity_matrix = np.exp(scaled_matrix)
    affinity_matrix[ix] = 0.  # mask
    affinity_matrix.flat[::len(affinity_matrix) + 1] = 0.  # diagonal
    return affinity_matrix


def double_centre(matrix, square_input=True):
    """ Double-centres the input matrix: From each element: Subtract the row
    mean Subtract the column mean Add the grand mean Divide by -2 Method
    from: Torgerson, W S (1952). Multidimensional scaling: I. Theory and
    method. Alternatively M = -0.5 * (I - 1/n)D[^2](I - 1/n) """

    m = matrix.copy()
    if square_input:
        m **= 2

    (rows, cols) = m.shape
    cm = np.mean(m, axis=0)  # column means
    rm = np.mean(m, axis=1).reshape((rows, 1))  # row means
    gm = np.mean(cm)  # grand mean
    m -= rm + cm - gm
    m /= -2
    return m

def _estimate_additive_constant(matrix):
    """
    CMDS Additive Constant: correction for non-Euclidean distances.
    Procedure taken from R function cmdscale.
    The additive constant is given by the largest eigenvalue (real part)
    of this 2x2 block matrix -

    /-------+-------\
    |       |       |
    |   0   | 2*dbc |    NB: dbc function(m: matrix):
    |       | (d^2) |        double_centre() [see above]
    +-------+-------+
    |       |       |
    |  -I   |-4*dbc |
    |       |  (2d) |
    \-------+-------/

    corrected matrix = matrix + additive constant (diagonal kept as 0)
    """
    topleft = np.zeros(matrix.shape)
    topright = 2*double_centre(matrix)
    bottomleft = -np.eye(matrix.shape[0])
    bottomright = -4*double_centre(matrix, square_input=False)
    Z = np.vstack([np.hstack([topleft,topright]),
                   np.hstack([bottomleft,bottomright])])
    return max(np.real(np.linalg.eigvals(Z)))


def _additive_correct(matrix):
    addc = _estimate_additive_constant(matrix)
    tmp = matrix + addc
    np.fill_diagonal(tmp, 0)
    return tmp


def check_euclidean(matrix):
    """ A distance matrix is euclidean iff F = -0.5 * (I - 1/n)D(I - 1/n) is
    PSD, where I is the identity matrix D is the distance matrix, `1` is a
    square matrix of ones, and n is the matrix size, common to all """

    f = double_centre(matrix, square_input=True)
    return check_psd(f)


def check_pd(matrix):
    """ A symmetric matrix (M) is PD if it has a Cholesky decomposition, i.e.
    M = R.T dot R, where R is upper triangular with positive diagonal entries
    """
    try:
        np.linalg.cholesky(matrix)
        return True
    except np.linalg.LinAlgError:
        return False


def check_psd(matrix, tolerance=1e-6):
    """ A square matrix is PSD if all eigenvalues of its Hermitian part are
    non- negative. The Hermitian part is given by (self + M*)/2, where M* is
    the complex conjugate transpose of M """

    hermitian = (matrix + matrix.T.conjugate()) / 2
    eigenvalues = np.linalg.eigh(hermitian)[0]
    return (eigenvalues > -tolerance).all()


def normalise_rows(matrix):
    """ Scales all rows to length 1. Fails when row is 0-length, so it
    leaves these unchanged """

    lengths = np.apply_along_axis(np.linalg.norm, 1, matrix)
    if not (lengths > 0).all():
        # raise ValueError('Cannot normalise 0 length vector to length 1')
        # print(matrix)
        lengths[lengths == 0] = 1
    return matrix / lengths[:, np.newaxis]


def kdists(matrix, k=7, ix=None):
    """ Returns the k-th nearest distances, row-wise, as a column vector """

    ix = ix or kindex(matrix, k)
    return matrix[ix][np.newaxis].T


def kindex(matrix, k):
    """ Returns indices to select the kth nearest neighbour"""

    ix = (np.arange(len(matrix)), matrix.argsort(axis=0)[k])
    return ix


def kmask(matrix, k=7, dists=None, logic='or'):
    """ Creates a boolean mask to include points within k nearest
    neighbours, and exclude the rest.
    Logic can be OR or AND. OR gives the k-nearest-neighbour mask,
    AND gives the mutual k-nearest-neighbour mask."""

    dists = (kdists(matrix, k=k) if dists is None else dists)
    mask = (matrix <= dists)
    if logic == 'or' or logic == '|':
        return mask | mask.T
    elif logic == 'and' or logic == '&':
        return mask & mask.T
    return mask


def kscale(matrix, k=7, dists=None):
    """ Returns the local scale based on the k-th nearest neighbour """
    dists = (kdists(matrix, k=k) if dists is None else dists)
    scale = dists.dot(dists.T)
    return scale


def laplace(affinity_matrix, shi_malik_type=False):
    """ Converts affinity matrix into normalised graph Laplacian,
    for spectral clustering.
    (At least) two forms exist:

    L = (D^-0.5).A.(D^-0.5) - default

    L = (D^-1).A - `Shi-Malik` type, from Shi Malik paper"""

    diagonal = affinity_matrix.sum(axis=1) - affinity_matrix.diagonal()
    zeros = diagonal <= 1e-10
    diagonal[zeros] = 1
    if (diagonal <= 1e-10).any():  # arbitrarily small value
        raise ZeroDivisionError
    if shi_malik_type:
        inv_d = np.diag(1 / diagonal)
        return inv_d.dot(affinity_matrix)
    diagonal = np.sqrt(diagonal)
    return affinity_matrix / diagonal / diagonal[:, np.newaxis]


def rbf(matrix, sigma=1):
    """ Returns the Radial Basis Function kernel from the distance matrix
    """
    return np.exp(-matrix ** 2 / (2 * sigma ** 2))


def shift_and_scale(matrix, shift, scale):
    """ Shift and scale matrix so its minimum value is placed at `shift` and
    its maximum value is scaled to `scale` """

    zeroed = matrix - matrix.min()
    scaled = (scale - shift) * (zeroed / zeroed.max())
    return scaled + shift


def get_permutation_matrix(input_ordering, desired_ordering):
    length = len(input_ordering)
    if not len(desired_ordering) == length:
        print('List lengths don\'t match')
        return
    p = np.zeros((length, length), dtype=np.int)
    for i in range(length):
        j = desired_ordering.index(input_ordering[i])
        p[i, j] = 1
    return p


def binsearch_dists(matrix, tolerance=0.001):
    mink = 1
    maxk = matrix.shape[0]
    guessk = int(np.log(maxk).round())
    last_result = (guessk, None)
    while maxk - mink != 1:
        dists = kdists(matrix, k=guessk)
        if (dists > tolerance).all():
            maxk = guessk
            last_result = (guessk, dists)
            guessk = mink + old_div((guessk - mink), 2)
        else:
            mink = guessk
            guessk += old_div((maxk - guessk), 2)

    if last_result[0] == guessk + 1:
        return last_result
    dists = kdists(matrix, k=guessk + 1)
    return guessk + 1, dists


def binsearch_mask(matrix, logic='or'):
    mink = 1
    maxk = matrix.shape[0]
    guessk = int(np.log(maxk).round())
    result = [guessk, None]
    while maxk - mink != 1:
        test_dist = kdists(matrix, k=guessk)
        test_mask = kmask(matrix, dists=test_dist, logic=logic)
        if isconnected(test_mask) and (test_dist > 1e-6).all():
            maxk = guessk  # either correct or too high
            result = (guessk, test_mask)
            guessk = mink + old_div((guessk - mink), 2)  # try a lower number
        else:
            mink = guessk  # too low
            guessk += old_div((maxk - guessk), 2)

    if result[0] == guessk + 1:
        k, mask = result
        scale = kscale(matrix, k)
        return k, mask, scale
    else:
        k = guessk + 1
        mask = kmask(matrix, k=guessk + 1, logic=logic)
        scale = kscale(matrix, k)
    return k, mask, scale


def eigen(matrix):
    """ Calculates the eigenvalues and eigenvectors of the input matrix.
    Returns a tuple of (eigenvalues, eigenvectors, cumulative percentage of
    variance explained). Eigenvalues and eigenvectors are sorted in order of
    eigenvalue magnitude, high to low """

    (vals, vecs) = np.linalg.eigh(matrix)
    ind = vals.argsort()[::-1]
    vals = vals[ind]
    vecs = vecs[:, ind]
    vals_ = vals.copy()
    vals_[vals_ < 0] = 0.
    cum_var_exp = np.cumsum(vals_ / vals_.sum())
    return Decomp(matrix.copy(), vals, vecs, cum_var_exp)


def _embedding_classical_mds(matrix, dimensions=3, additive_correct=False):
    """
    Private method to calculate CMDS embedding
    :param dimensions: (int)
    :return: coordinate matrix (np.array)
    """
    if additive_correct:
        dbc = double_centre(_additive_correct(matrix))
    else:
        dbc = double_centre(matrix)
    decomp = eigen(dbc)
    lambda_ = np.diag(np.sqrt(np.abs(decomp.vals[:dimensions])))
    evecs = decomp.vecs[:, :dimensions]
    coords = evecs.dot(lambda_)
    return coords


def _embedding_spectral(matrix, dimensions=3, unit_length=True,
                        affinity_matrix=None, sigma=1):
    """
    Private method to calculate Spectral embedding
    :param dimensions: (int)
    :return: coordinate matrix (np.array)
    """
    if affinity_matrix is None:
        aff = rbf(matrix, sigma=sigma)
    else:
        aff = affinity_matrix
    coords = sklearn.manifold.spectral_embedding(aff, dimensions)
    return normalise_rows(coords) if unit_length else coords


def _embedding_tsne(matrix, dimensions=3, early_exaggeration=12.0,
                    method='barnes_hut', perplexity=30, learning_rate=200,
                    n_iter=1000):
    """
    Private method to perform tSNE embedding
    :param matrix: treeCl Distance Matrix
    :param dimensions: Number of dimensions in which to embed points
    :return: treeCl CoordinateMatrix
    """
    tsne = sklearn.manifold.TSNE(n_components=dimensions,
                                 metric="precomputed",
                                 early_exaggeration=early_exaggeration,
                                 method=method,
                                 perplexity=perplexity,
                                 learning_rate=learning_rate,
                                 n_iter=1000)
    return tsne.fit_transform(matrix)


def _embedding_metric_mds(matrix, dimensions=3):
    """
    Private method to calculate MMDS embedding
    :param dimensions: (int)
    :return: coordinate matrix (np.array)
    """
    mds = sklearn.manifold.MDS(n_components=dimensions,
                               dissimilarity='precomputed',
                               metric=True)
    mds.fit(matrix)
    return mds.embedding_


def _embedding_nonmetric_mds(matrix, dimensions=3, initial_coords=None):
    """
    Private method to calculate NMMDS embedding
    :param dimensions: (int)
    :return: coordinate matrix (np.array)
    """
    mds = sklearn.manifold.MDS(n_components=dimensions,
                               dissimilarity='precomputed',
                               metric=False)
    if initial_coords is not None:
        mds.fit(matrix, init=initial_coords)
    else:
        mds.fit(matrix)
    return mds.embedding_


def _embedding_kernel_pca(matrix, dimensions=3, affinity_matrix=None,
                          sigma=1):
    """
    Private method to calculate KPCA embedding
    :param dimensions: (int)
    :return: coordinate matrix (np.array)
    """
    if affinity_matrix is None:
        aff = rbf(matrix, sigma)
    else:
        aff = affinity_matrix
    kpca = sklearn.decomposition.KernelPCA(kernel='precomputed',
                                           n_components=dimensions)
    return kpca.fit_transform(aff)


class Decomp(object):
    """ Eigen decomposition result """

    def __init__(self, matrix, vals, vecs, cve):
        self.matrix = matrix
        self.vals = vals
        self.vecs = vecs
        self.cve = cve

    def __str__(self):
        return '\n'.join([str(self.vals), str(self.vecs), str(self.cve)])

    def coords_by_cutoff(self, cutoff=0.80):
        """ Returns fitted coordinates in as many dimensions as are needed to
        explain a given amount of variance (specified in the cutoff) """

        i = np.where(self.cve >= cutoff)[0][0]
        coords_matrix = self.vecs[:, :i + 1]
        return coords_matrix, self.cve[i]

    def coords_by_dimension(self, dimensions=3):
        """ Returns fitted coordinates in specified number of dimensions, and
        the amount of variance explained) """

        coords_matrix = self.vecs[:, :dimensions]
        varexp = self.cve[dimensions - 1]
        return coords_matrix, varexp


class Matrix(object):
    def __repr__(self):
        return repr(self.df)

    @property
    def values(self):
        return self.to_array()

    @property
    def shape(self):
        return self.to_array().shape

    def to_array(self):
        return self.df.values


class CoordinateMatrix(Matrix):
    def __init__(self, array, names=None):
        nrow = array.shape[0]
        if names is not None and len(names) != nrow:
            msg = 'Length of names array doesn\'t match the number of rows in the coordinate matrix ({}, {})'.format(
                len(names), nrow)
            raise ValueError(msg)

        self.df = pd.DataFrame(array, index=names)


class DistanceMatrix(Matrix):
    def __init__(self):
        self.df = pd.DataFrame()

    def __eq__(self, other):
        if (np.abs(self.sort().values - other.sort().values) < 1e-10).all():
            return True

    @classmethod
    def from_csv(cls, filename, **kwargs):
        with fileIO.freader(filename) as handle:
            new_instance = cls()
            new_instance.df = pd.DataFrame.from_csv(handle, **kwargs)
        return new_instance

    def to_csv(self, filename, **kwargs):
        with fileIO.fwriter(filename) as handle:
            self.df.to_csv(handle)

    @classmethod
    def from_array(cls, array, names=None):
        nrow, ncol = array.shape
        if nrow != ncol:
            errmsg = 'Not a square array'
            raise ValueError(errmsg)
        new_instance = cls()
        new_instance.df = pd.DataFrame(array)
        try:
            new_instance.set_names(names)
            return new_instance
        except ValueError as err:
            sys.stderr.write(str(err))
            return new_instance

    def get_names(self):
        return [str(x) for x in self.df.index]

    def set_names(self, names):
        if names is None:
            return
        if len(names) != self.df.index.size:
            errmsg = 'Expected {} names, got {}'.format(self.df.index.size, len(names))
            raise ValueError(errmsg)
        self.df.index = self.df.columns = names

    def add_noise(self, std=0.0001):
        ix = np.triu_indices(len(self.df), 1)
        rev_ix = ix[::-1]
        noise = np.random.normal(0, std, len(ix[0]))
        noisy = self.to_array().copy()
        noisy[ix] += noise
        noisy[rev_ix] += noise
        noisy[noisy < 0] = np.abs(noisy[noisy < 0])
        return self.__class__(noisy, self.df.index)

    def affinity(self, mask=None, scale=None):
        return affinity(self.to_array(), mask, scale)

    def double_centre(self, square_input=True):
        return double_centre(self.to_array(), square_input)

    def embedding(self, dimensions, method, **kwargs):
        """
        Embeds the distance matrix in a coordinate space. Implemented methods are:
            cmds: Classical MultiDimensional Scaling
            kpca: Kernel Principal Components Analysis
            mmds: Metric MultiDimensional Scaling
            nmmds: Non-Metric MultiDimensional Scaling
            spectral: Spectral decomposition of Laplacian matrix
            tsne: t-distributed Stochastic Neighbour Embedding

        Valid kwargs:
            kpca: affinity_matrix - a precomputed array of affinities
                  sigma - the value of sigma to use when computing the affinity matrix via
                          the Radial Basis Function
            nmmds: initial_coords - a set of coordinates to refine. NMMDS works very badly
                                    without this
            spectral: affinity_matrix, sigma
                      unit_length - scale the coordinates to unit length, so points sit
                                    on the surface of the unit sphere
        :param dimensions: (int) number of coordinate axes to use
        :param method: (string) one of cmds, kpca, mmds, nmmds, spectral
        :param kwargs: unit_length (bool), affinity_matrix (np.array), sigma (float), initial_coords (np.array)
        :return: coordinate matrix (np.array)
        """
        errors.optioncheck(method, ['cmds', 'kpca', 'mmds', 'nmmds', 'spectral', 'tsne'])
        if method == 'cmds':
            array =  _embedding_classical_mds(self.to_array(), dimensions, **kwargs)
        elif method == 'kpca':
            array = _embedding_kernel_pca(self.to_array(), dimensions, **kwargs)
        elif method == 'mmds':
            array = _embedding_metric_mds(self.to_array(), dimensions)
        elif method == 'nmmds':
            array = _embedding_nonmetric_mds(self.to_array(), dimensions, **kwargs)
        elif method == 'spectral':
            array = _embedding_spectral(self.to_array(), dimensions, **kwargs)
        elif method == 'tsne':
            array = _embedding_tsne(self.to_array(), dimensions, **kwargs)

        return CoordinateMatrix(array, names=self.df.index)

    def reorder(self, new_order):
        reordered_df = self.df.reindex(columns=new_order, index=new_order)
        reordered_names = reordered_df.columns
        newobj = self.__class__()
        newobj.df = reordered_df
        return newobj

    def sort(self):
        order = self.df.index.argsort()
        return self.reorder(order)
