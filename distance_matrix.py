#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm as CM
from externals import GTP
import utils


def get_dendropy_distances(dpy_trees, fn):
    num_trees = len(dpy_trees)
    matrix = np.zeros((num_trees, num_trees))
    for i in range(num_trees):
        for j in range(i + 1, num_trees):
            distance = fn(dpy_trees[i], dpy_trees[j])
            matrix[i, j] = matrix[j, i] = distance
    return matrix


def get_geo_distances(trees, tmpdir=None):

    g = GTP(tmpdir=tmpdir)
    return g.run(trees)


def get_distance_matrix(trees, metric, tmpdir):
    """ Generates pairwise distance matrix between trees Uses one of the
    following distance metrics: Robinson-Foulds distance - topology only (='rf')
    Robinson-Foulds distance - branch lengths (='wrf') Euclidean distance -
    Felsenstein's branch lengths distance (='euc') Geodesic distance - branch
    lengths (='geo') """

    if metric == 'geo':
        return get_geo_distances(trees, tmpdir=tmpdir)

    dpy_trees = utils.dpy.convert_to_dendropy_trees(trees)
    if metric == 'rf':
        matrix = get_dendropy_distances(dpy_trees, utils.dpy.get_rf_distance)
    elif metric == 'wrf':

        matrix = get_dendropy_distances(dpy_trees, utils.dpy.get_wrf_distance)
    elif metric == 'euc':

        matrix = get_dendropy_distances(dpy_trees, utils.dpy.get_euc_distance)
    else:

        print 'Unrecognised distance metric'
        return
    return matrix


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


class Decomp(object):

    """ Eigen decomposition result """

    def __init__(
        self,
        matrix,
        vals,
        vecs,
        cve,
        ):
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
        return (coords_matrix, self.cve[i])

    def coords_by_dimension(self, dimensions=3):
        """ Returns fitted coordinates in specified number of dimensions, and
        the amount of variance explained) """

        coords_matrix = self.vecs[:, :dimensions]
        varexp = self.cve[dimensions - 1]
        return (coords_matrix, varexp)


class DistanceMatrix(np.ndarray):

    def __new__(
        cls,
        trees,
        metric,
        tmpdir='/tmp',
        dtype=float,
        ):

        input_array = get_distance_matrix(trees, metric, tmpdir)
        obj = np.asarray(input_array, dtype).view(cls)
        obj.metric = metric
        obj.tmpdir = tmpdir
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.metric = getattr(obj, 'metric', None)
        self.tmpdir = getattr(obj, 'tmpdir', None)

    def __array_wrap__(self, out_arr, context=None):

        # print context

        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __eq__(self, other):
        if (np.abs(self - other) < 1e-10).all():
            return True

    def add_noise(self):
        ix = np.triu_indices(len(self), 1)
        rev_ix = ix[::-1]
        noise = np.random.normal(0, 0.001, len(ix[0]))
        noisy = self.copy()
        noisy[ix] += noise
        noisy[rev_ix] += noise
        noisy[noisy < 0] = np.abs(noisy[noisy < 0])
        return noisy

    def affinity(
        self,
        mask=None,
        scale=None,
        sigma=2,
        ):

        mask = (mask if mask is not None else np.ones(self.shape, dtype=bool))
        scale = (scale if scale is not None else np.ones(self.shape) * 2
                 * sigma)
        ix = np.where(np.logical_not(mask))
        affinity_matrix = np.exp(-self ** 2 / scale)
        affinity_matrix[ix] = 0.  # mask
        affinity_matrix.flat[::len(affinity_matrix) + 1] = 0.  # diagonal
        return affinity_matrix

    def binsearch_dists(self):
        mink = 1
        maxk = self.shape[0]
        guessk = int(np.log(maxk).round())
        last_result = (guessk, None)
        while maxk - mink != 1:
            dists = self.kdists(k=guessk)
            if (dists > 0).all():
                maxk = guessk
                last_result = (guessk, dists)
                guessk = mink + (guessk - mink) / 2
            else:
                mink = guessk
                guessk = guessk + (maxk - guessk) / 2

        if last_result[0] == guessk + 1:
            return last_result
        dists = self.kdists(k=guessk + 1)
        return (guessk + 1, dists)

    def binsearch_mask(self, logic='or'):

        mink = 1
        maxk = self.shape[0]
        guessk = int(np.log(maxk).round())
        result = [guessk, None]
        while maxk - mink != 1:
            test = self.kmask(k=guessk, logic=logic)
            if isconnected(test):
                maxk = guessk  # either correct or too high
                result = (guessk, test)
                guessk = mink + (guessk - mink) / 2  # try a lower number
            else:
                mink = guessk  # too low
                guessk = guessk + (maxk - guessk) / 2

        if result[0] == guessk + 1:
            return result
        mask = self.kmask(k=guessk + 1, logic=logic)
        return (guessk + 1, mask)

    def check_euclidean(self):
        """ A distance matrix is euclidean iff F = -0.5 * (I - 1/n)D(I - 1/n) is
        PSD, where I is the identity matrix D is the distance matrix, `1` is a
        square matrix of ones, and n is the matrix size, common to all """

        F = self.double_centre(square_input=True)
        return F.check_PSD()

    def check_PD(self):
        """ A symmetric matrix is PD if 1: all diagonal entries are positive,
        and 2: the diagonal element is greater than the sum of all other entries
        """

        diagonal = self.diagonal()
        colsum = self.sum(axis=0)
        rowsum = self.sum(axis=1)
        symmetric = self.T == self  # Predicates
        diag_pos = (diagonal > 0).all()
        col_dominant = (diagonal > colsum - diagonal).all()
        row_dominant = (diagonal > rowsum - diagonal).all()
        return np.all([symmetric, diag_pos, col_dominant, row_dominant])

    def check_PSD(self):
        """ A square matrix is PSD if all eigenvalues of its Hermitian part are
        non- negative. The Hermitian part is given by (self + M*)/2, where M* is
        the complex conjugate transpose of M """

        hermitian = (self + self.T.conjugate()) / 2
        eigenvalues = np.linalg.eigh(hermitian)[0]
        return not (eigenvalues < 0).all()

    def double_centre(self, square_input=True):
        """ Double-centres the input matrix: From each element: Subtract the row
        mean Subtract the column mean Add the grand mean Divide by -2 Method
        from: Torgerson, W S (1952). Multidimensional scaling: I. Theory and
        method. Alternatively F = -0.5 * (I - 1/n)D[^2](I - 1/n) """

        M = (self * self if square_input else self.copy())
        (rows, cols) = M.shape

        cm = np.mean(M, axis=0)  # column means
        rm = np.mean(M, axis=1).reshape((rows, 1))  # row means
        gm = np.mean(cm)  # grand mean
        M -= rm + cm - gm
        M /= -2
        return M

    def eigen(self):
        """ Calculates the eigenvalues and eigenvectors of the input matrix.
        Returns a tuple of (eigenvalues, eigenvectors, cumulative percentage of
        variance explained). Eigenvalues and eigenvectors are sorted in order of
        eigenvalue magnitude, high to low """

        (vals, vecs) = np.linalg.eigh(self)
        ind = vals.argsort()[::-1]
        vals = vals[ind]
        vecs = vecs[:, ind]
        vals_ = vals.copy()
        vals_[vals_ < 0] = 0.
        cum_var_exp = np.cumsum(vals_ / vals_.sum())
        return Decomp(self.copy, vals, vecs, cum_var_exp)

    def normalise_rows(self):
        """ Scales all rows to length 1. Fails when row is 0-length, so it
        leaves these unchanged """

        lengths = np.apply_along_axis(np.linalg.norm, 1, self)
        if not (lengths > 0).all():
            print 'Cannot normalise 0 length vector to length 1'
            print self
            lengths[lengths == 0] = 1
        return self / lengths[:, np.newaxis]

    def kdists(self, k=7, ix=None):
        """ Returns the k-th nearest distances """

        ix = ix or self.kindex(k)
        return self[ix]

    def kindex(self, k):
        """ Returns indices to select the kth nearest neighbour, column-wise"""

        ix = list(np.ix_(*[np.arange(i) for i in self.shape]))
        ix[0] = self.argsort(0)[k - 1:k, :]
        return ix

    def kmask(
        self,
        dists=None,
        k=7,
        logic='or',
        ):
        """ Creates a boolean mask to include points within k nearest
        neighbours, and exclude the rest """

        dists = (self.kdists(k=k) if dists is None else dists)
        if logic == 'or' or logic == '|':
            mask = (self <= dists) | (self <= dists.T)
        elif logic == 'and' or logic == '&':
            mask = (self <= dists) & (self <= dists.T)
        return mask

    def kscale(self, dists=None, k=7):
        """ Returns the local scale based on the k-th nearest neighbour """

        dists = (self.kdists(k=k) if dists is None else dists)
        scale = dists.T.dot(dists)
        return scale

    def laplace(self, affinity_matrix, shi_malik_type=False):
        """ Converts affinity matrix into graph Laplacian, for spectral
        clustering. (At least) two forms exist:
        
        L = (D^-0.5).A.(D^-0.5) - default
        
        L = (D^-1).A - `Shi-Malik` type, from Shi Malik paper"""

        diagonal = affinity_matrix.sum(axis=1) - affinity_matrix.diagonal()
        if 0. in diagonal:
            raise ZeroDivisionError
        if shi_malik_type:
            invD = np.diag(1 / diagonal).view(type(self))
            return invD.dot(affinity_matrix)
        invRootD = np.diag(np.sqrt(1 / diagonal)).view(type(self))
        return invRootD.dot(affinity_matrix).dot(invRootD)

    def normalise(self):
        """ Shift and scale matrix to [0,1] interval """

        return self.shift_and_scale(0, 1)

    def shift_and_scale(self, shift, scale):
        """ Shift and scale matrix so its minimum value is placed at `shift` and
        its maximum value is scaled to `scale` """

        zeroed = self - self.min()
        scaled = scale * (zeroed / zeroed.max())
        return scaled + shift

    def get_permutation_matrix(self, input_ordering, desired_ordering):
        length = len(input_ordering)
        if not len(desired_ordering) == length:
            print 'List lengths don\'t match'
            return
        P = np.zeros((length, length), dtype=np.int)
        for i in range(length):
            j = desired_ordering.index(input_ordering[i])
            P[i, j] = 1
        return P

    def plot_heatmap(self, sort_partition=None):
        """ Sort partition should be a flatlist of the clusters as returned by
        Partition().get_memberships(..., flatten=True) """

        length = self.shape[0]
        datamax = float(np.abs(self).max())
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ticks_at = [0, 0.5 * datamax, datamax]
        if sort_partition:
            p = self.get_permutation_matrix(range(len(sort_partition)),
                    sort_partition)
            self = np.dot(p.T, np.dot(self, p))
        cax = ax.imshow(
            np.array(self),
            interpolation='nearest',
            origin='lower',
            extent=[0., length, 0., length],
            vmin=0,
            vmax=datamax,
            cmap=CM.Blues,
            )
        cbar = fig.colorbar(cax, ticks=ticks_at, format='%1.2g')
        cbar.set_label('Distance')
        return fig
