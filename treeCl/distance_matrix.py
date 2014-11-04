#!/usr/bin/env python
from __future__ import print_function

# third party
import numpy as np
import sklearn

# treeCl
import errors
from treedist import eucdist_matrix, geodist_matrix, rfdist_matrix, wrfdist_matrix

dist_mtx_fns = dict(euc=eucdist_matrix, geo=geodist_matrix, rf=rfdist_matrix, wrf=wrfdist_matrix)


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
        return coords_matrix, self.cve[i]

    def coords_by_dimension(self, dimensions=3):
        """ Returns fitted coordinates in specified number of dimensions, and
        the amount of variance explained) """

        coords_matrix = self.vecs[:, :dimensions]
        varexp = self.cve[dimensions - 1]
        return coords_matrix, varexp


# noinspection PyNoneFunctionAssignment
class DistanceMatrix(np.ndarray):
    # noinspection PyNoneFunctionAssignment
    def __new__(
            cls,
            trees,
            metric,
            dtype=float,
            add_noise=False,
            normalise=False,
    ):
        errors.optioncheck(metric, ['euc', 'geo', 'rf', 'wrf'])
        fn = dist_mtx_fns[metric]
        input_array = fn(trees, normalise)
        obj = np.asarray(input_array, dtype).view(cls)
        obj.metric = metric
        if add_noise:
            obj = obj.add_noise()
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.metric = getattr(obj, 'metric', None)
        self.tmpdir = getattr(obj, 'tmpdir', None)

    def __array_wrap__(self, out_arr, context=None):

        # print(context)

        # noinspection PyArgumentList
        return np.ndarray.__array_wrap__(self, out_arr, context)

    def __eq__(self, other):
        if (np.abs(self - other) < 1e-10).all():
            return True

    def add_noise(self):
        ix = np.triu_indices(len(self), 1)
        rev_ix = ix[::-1]
        noise = np.random.normal(0, 0.0001, len(ix[0]))
        noisy = self.copy()
        noisy[ix] += noise
        noisy[rev_ix] += noise
        noisy[noisy < 0] = np.abs(noisy[noisy < 0])
        return noisy

    def affinity(
            self,
            mask=None,
            scale=None,
    ):
        """
        Mask is a 2d boolean matrix. Scale is a 2d local scale matrix,
        as output by self.kscale(). It's the outer product of the kdists column
        vector produced by self.kdists.
        """

        mask = (mask if mask is not None else np.ones(self.shape, dtype=bool))
        assert isconnected(mask)
        scale = (scale if scale is not None else np.ones(self.shape))

        ix = np.where(np.logical_not(mask))
        scaled_matrix = -self ** 2 / scale
        # inputs where distance = 0 and scale = 0 result in NaN:
        # the next line replaces NaNs with -1.0
        scaled_matrix[np.where(np.isnan(scaled_matrix))] = -1.0
        affinity_matrix = np.exp(scaled_matrix)
        affinity_matrix[ix] = 0.  # mask
        affinity_matrix.flat[::len(affinity_matrix) + 1] = 0.  # diagonal
        return affinity_matrix

    def binsearch_dists(self, tolerance=0.001):
        mink = 1
        maxk = self.shape[0]
        guessk = int(np.log(maxk).round())
        last_result = (guessk, None)
        while maxk - mink != 1:
            dists = self.kdists(k=guessk)
            if (dists > tolerance).all():
                maxk = guessk
                last_result = (guessk, dists)
                guessk = mink + (guessk - mink) / 2
            else:
                mink = guessk
                guessk += (maxk - guessk) / 2

        if last_result[0] == guessk + 1:
            return last_result
        dists = self.kdists(k=guessk + 1)
        return guessk + 1, dists

    def binsearch_mask(self, logic='or'):

        mink = 1
        maxk = self.shape[0]
        guessk = int(np.log(maxk).round())
        result = [guessk, None]
        while maxk - mink != 1:
            test_dist = self.kdists(k=guessk)
            test_mask = self.kmask(dists=test_dist, logic=logic)
            if isconnected(test_mask) and (test_dist > 1e-6).all():
                maxk = guessk  # either correct or too high
                result = (guessk, test_mask)
                guessk = mink + (guessk - mink) / 2  # try a lower number
            else:
                mink = guessk  # too low
                guessk += (maxk - guessk) / 2

        if result[0] == guessk + 1:
            k, mask = result
            scale = self.kscale(k)
            return k, mask, scale
        else:
            k = guessk + 1
            mask = self.kmask(k=guessk + 1, logic=logic)
            scale = self.kscale(k)
        return k, mask, scale

    def check_euclidean(self):
        """ A distance matrix is euclidean iff F = -0.5 * (I - 1/n)D(I - 1/n) is
        PSD, where I is the identity matrix D is the distance matrix, `1` is a
        square matrix of ones, and n is the matrix size, common to all """

        f = self.double_centre(square_input=True)
        return f.check_psd()

    def check_pd(self):
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

    def check_psd(self):
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
        method. Alternatively M = -0.5 * (I - 1/n)D[^2](I - 1/n) """

        m = (self * self if square_input else self.copy())
        (rows, cols) = m.shape

        cm = np.mean(m, axis=0)  # column means
        rm = np.mean(m, axis=1).reshape((rows, 1))  # row means
        gm = np.mean(cm)  # grand mean
        m -= rm + cm - gm
        m /= -2
        return m

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
        return Decomp(self.copy(), vals, vecs, cum_var_exp)

    def embedding(self, dimensions, method, **kwargs):
        """
        Embeds the distance matrix in a coordinate space. Implemented methods are:
            cmds: Classical MultiDimensional Scaling
            kpca: Kernel Principal Components Analysis
            mmds: Metric MultiDimensional Scaling
            nmmds: Non-Metric MultiDimensional Scaling
            spectral: Spectral decomposition of Laplacian matrix

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
        errors.optioncheck(method, ['cmds', 'kpca', 'mmds', 'nmmds', 'spectral'])
        if method == 'cmds':
            return self._embedding_classical_mds(dimensions)
        elif method == 'kpca':
            return self._embedding_kernel_pca(dimensions, **kwargs)
        elif method == 'mmds':
            return self._embedding_metric_mds(dimensions)
        elif method == 'nmmds':
            return self._embedding_nonmetric_mds(dimensions, **kwargs)
        elif method == 'spectral':
            return self._embedding_spectral(dimensions, **kwargs)

    def _embedding_classical_mds(self, dimensions=3):
        """
        Private method to calculate CMDS embedding
        :param dimensions: (int)
        :return: coordinate matrix (np.array)
        """
        dbc = self.double_centre()
        decomp = dbc.eigen()
        lambda_ = np.diag(np.sqrt(np.abs(decomp.vals[:dimensions])))
        evecs = decomp.vecs[:, :dimensions]
        coords = evecs.dot(lambda_)
        return coords

    def _embedding_spectral(self, dimensions=3, unit_length=True,
                            affinity_matrix=None, sigma=1):
        """
        Private method to calculate Spectral embedding
        :param dimensions: (int)
        :return: coordinate matrix (np.array)
        """
        if affinity_matrix is None:
            aff = self.rbf(sigma=sigma)
        else:
            aff = affinity_matrix
        coords = aff.laplace(aff).eigen().vecs[:, :dimensions]
        if unit_length:  # normalise all vectors to unit length
            coords /= np.sqrt((coords ** 2).sum(axis=1))[:, np.newaxis]
        return coords

    def _embedding_metric_mds(self, dimensions=3):
        """
        Private method to calculate MMDS embedding
        :param dimensions: (int)
        :return: coordinate matrix (np.array)
        """
        mds = sklearn.manifold.MDS(n_components=dimensions,
                                   dissimilarity='precomputed',
                                   metric=True)
        mds.fit(self)
        return mds.embedding_

    def _embedding_nonmetric_mds(self, dimensions=3, initial_coords=None):
        """
        Private method to calculate NMMDS embedding
        :param dimensions: (int)
        :return: coordinate matrix (np.array)
        """
        mds = sklearn.manifold.MDS(n_components=dimensions,
                                   dissimilarity='precomputed',
                                   metric=False)
        if initial_coords is not None:
            mds.fit(self, init=initial_coords)
        else:
            mds.fit(self)
        return mds.embedding_

    def _embedding_kernel_pca(self, dimensions=3, affinity_matrix=None,
                              sigma=1):
        """
        Private method to calculate KPCA embedding
        :param dimensions: (int)
        :return: coordinate matrix (np.array)
        """
        if affinity_matrix is None:
            aff = self.rbf(sigma)
        else:
            aff = affinity_matrix
        kpca = sklearn.decomposition.KernelPCA(kernel='precomputed',
                                               n_components=dimensions)
        return kpca.fit_transform(aff)

    def normalise_rows(self):
        """ Scales all rows to length 1. Fails when row is 0-length, so it
        leaves these unchanged """

        lengths = np.apply_along_axis(np.linalg.norm, 1, self)
        if not (lengths > 0).all():
            print('Cannot normalise 0 length vector to length 1')
            print(self)
            lengths[lengths == 0] = 1
        return self / lengths[:, np.newaxis]

    def kdists(self, k=7, ix=None):
        """ Returns the k-th nearest distances, row-wise, as a column vector """

        ix = ix or self.kindex(k)
        return self[ix][np.newaxis].T

    def kindex(self, k):
        """ Returns indices to select the kth nearest neighbour"""

        ix = (np.arange(len(self)), self.argsort(axis=0)[k])
        # ix = list(np.ix_(*[np.arange(i) for i in self.shape]))
        # ix[0] = self.argsort(0)[k:k+1, :]
        return ix

    def kmask(
            self,
            k=7,
            dists=None,
            logic='or',
    ):
        """ Creates a boolean mask to include points within k nearest
        neighbours, and exclude the rest.
        Logic can be OR or AND. OR gives the k-nearest-neighbour mask,
        AND gives the mutual k-nearest-neighbour mask."""

        dists = (self.kdists(k=k) if dists is None else dists)
        mask = (self <= dists)
        if logic == 'or' or logic == '|':
            return mask | mask.T
        elif logic == 'and' or logic == '&':
            return mask & mask.T
        return mask

    def kscale(self, k=7, dists=None):
        """ Returns the local scale based on the k-th nearest neighbour """

        dists = (self.kdists(k=k) if dists is None else dists)
        scale = dists.dot(dists.T)
        return scale

    def laplace(self, affinity_matrix, shi_malik_type=False):
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
            inv_d = np.diag(1 / diagonal).view(type(self))
            return inv_d.dot(affinity_matrix)
        diagonal = np.sqrt(diagonal)
        return affinity_matrix / diagonal / diagonal[:, np.newaxis]

    def normalise(self):
        """ Shift and scale matrix to [0,1] interval """

        return self.shift_and_scale(0, 1)

    def rbf(self, sigma=1):
        """ Returns the Radial Basis Function kernel from the distance matrix
        """
        return np.exp(-self ** 2 / (2 * sigma ** 2))

    def shift_and_scale(self, shift, scale):
        """ Shift and scale matrix so its minimum value is placed at `shift` and
        its maximum value is scaled to `scale` """

        zeroed = self - self.min()
        scaled = (scale - shift) * (zeroed / zeroed.max())
        return scaled + shift

    @staticmethod
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
