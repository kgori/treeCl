#!/usr/bin/env python
from __future__ import print_function

# standard library

# third party
import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from fastcluster import linkage
import skbio

try:
    from Bio.Cluster import kmedoids
    Biopython_Unavailable = False
except ImportError:
    print("Biopython unavailable - kmedoids clustering disabled")
    Biopython_Unavailable = True

from sklearn.cluster import AffinityPropagation, DBSCAN, KMeans
from sklearn.mixture import GMM

# treeCl
from distance_matrix import DistanceMatrix, rbf, binsearch_mask, kmask, kscale, affinity, laplace, eigen, double_centre, \
    normalise_rows
from partition import Partition


class Clustering(object):
    """ Apply clustering methods to distance matrix

    = Hierarchical clustering - single-linkage - complete-linkage - average-
    linkage (UPGMA) - Ward's method

    = k-medoids

    = Multidimensional Scaling (Principal Coordinate Analysis) + k-means

    = Spectral Clustering + k-means - NJW method - Shi-Malik method - Zelnik-
    Manor and Perona Local Scaling - Local Scaling with eigenvector rotation as
    stopping criterion

    """

    def __init__(self, dm):
        if isinstance(dm, np.ndarray):
            dm = DistanceMatrix(dm)

        if not isinstance(dm, DistanceMatrix):
            raise ValueError('Distance matrix should be a numpy array or treeCl.DistanceMatrix')
        self.dm = dm

    def __str__(self):
        return str(self.dm)

    def anosim(self, partition, n_permutations=999):
        if partition.is_minimal():
            raise ValueError("ANOSim is not defined for singleton clusters")
        elif partition.is_maximal():
            raise ValueError("ANOSim is not defined for maximally divided partitions")
        result = skbio.stats.distance.ANOSIM(skbio.DistanceMatrix(self.get_dm(False)), partition.partition_vector)
        return result(n_permutations)

    def permanova(self, partition, n_permutations=999):
        if partition.is_minimal():
            raise ValueError("PERMANOVA is not defined for singleton clusters")
        elif partition.is_maximal():
            raise ValueError("PERMANOVA is not defined for maximally divided partitions")
        result = skbio.stats.distance.PERMANOVA(skbio.DistanceMatrix(self.get_dm(False)), partition.partition_vector)
        return result(n_permutations)

    def kmedoids(self, nclusters, noise=False, npass=100, nreps=1):

        if Biopython_Unavailable:
            print('kmedoids not available without Biopython')
            return

        matrix = self.get_dm(noise)

        p = [kmedoids(matrix, nclusters=nclusters, npass=npass) for _ in
             range(nreps)]
        p.sort(key=lambda x: x[1])
        return Partition(p[0][0])

    def affinity_propagation(self, affinity_matrix=None, sigma=1, **kwargs):
        """

        :param kwargs: damping=0.5, max_iter=200, convergence_iter=15, copy=True, preference=None, verbose=False
        :return:
        """
        if affinity_matrix is None:
            aff = rbf(self.dm.values, sigma)
        else:
            aff = affinity_matrix

        est = AffinityPropagation(affinity='precomputed', **kwargs)
        est.fit(aff.view(np.ndarray))
        return Partition(est.labels_)

    def dbscan(self, eps=0.75, min_samples=3):
        """
        :param kwargs: key-value arguments to pass to DBSCAN
                       (eps: max dist between points in same neighbourhood,
                        min_samples: number of points in a neighbourhood)
        :return:
        """
        est = DBSCAN(metric='precomputed', eps=eps, min_samples=min_samples)
        est.fit(self.get_dm(False))
        return Partition(est.labels_)

    def get_dm(self, noise):
        return self.dm.add_noise().values if noise else self.dm.values

    def hierarchical(self, nclusters, linkage_method, noise=False):
        """

        :param nclusters: Number of clusters to return
        :param linkage_method: single, complete, average, ward, weighted, centroid or median
                               (http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html)
        :param noise: Add Gaussian noise to the distance matrix prior to clustering (bool, default=False)
        :return: Partition object describing clustering
        """
        matrix = self.get_dm(noise)

        linkmat = linkage(squareform(matrix), linkage_method)
        linkmat_size = len(linkmat)
        if nclusters <= 1:
            br_top = linkmat[linkmat_size - nclusters][2]
        else:
            br_top = linkmat[linkmat_size - nclusters + 1][2]
        if nclusters >= len(linkmat):
            br_bottom = 0
        else:
            br_bottom = linkmat[linkmat_size - nclusters][2]
        threshold = 0.5 * (br_top + br_bottom)
        t = fcluster(linkmat, threshold, criterion='distance')
        return Partition(t)

    def spectral_decomp(
            self,
            prune=None,
            local_scale=None,
            noise=False,
            verbosity=0,
            logic='or',
            **kwargs):
        """ Use prune to remove links between distant points:

        prune is None: no pruning
        prune={int > 0}: prunes links beyond `prune` nearest neighbours
        prune='estimate': searches for the smallest value that retains a fully
        connected graph

        """

        matrix = self.get_dm(noise)

        kp, mask, est_scale = binsearch_mask(matrix, logic=logic)  # prune anyway,

        # get local scale estimate
        ks = kp  # ks and kp are the scaling and pruning parameters
        est_ks = ks

        # ADJUST MASK
        if prune is None:  # change mask to all
            kp = len(matrix) - 1
            mask = np.ones(matrix.shape, dtype=bool)
        elif isinstance(prune, int) and prune > 0:
            kp = prune
            mask = kmask(matrix, prune, logic=logic)
        else:
            if not prune=='estimate':
                raise ValueError("'prune' should be None, a positive integer value, or 'estimate', not {}".format(prune))

        # ADJUST SCALE
        if local_scale is not None:
            if local_scale == 'median':
                ks = 'median'
                dist = np.median(matrix, axis=1)
                scale = np.outer(dist, dist)
            elif isinstance(local_scale, int):
                ks = local_scale
                scale = kscale(matrix, local_scale)
            else:
                scale = est_scale
        else:
            scale = est_scale

        # ZeroDivisionError safety check
        if not (scale > 1e-5).all():
            if verbosity > 0:
                print('Rescaling to avoid zero-div error')
            scale = est_scale
            ks = est_ks
            assert (scale > 1e-5).all()

        aff = affinity(matrix, mask, scale)
        laplacian = laplace(aff, **kwargs)

        if verbosity > 0:
            print('Pruning parameter: {0}'.format(kp))
            print('Scaling parameter: {0}'.format(ks))
            print('Mask, scale, affinity matrix and laplacian:')
            print(mask)
            print(scale)
            print(aff)
            print(laplace)

        return eigen(laplacian)  # vectors are in columns

    def spectral_cluster(self, nclusters, decomp, verbosity=0):

        if nclusters == 1:
            return Partition([1] * len(self.dm))

        pos = 0
        for val in decomp.vals:
            if val > 0:
                pos += 1

        (coords, cve) = decomp.coords_by_dimension(min(max(nclusters, 3), pos))
        if verbosity > 0:
            print('{0} dimensions explain {1:.2f}% of '
                  'the variance'.format(nclusters, cve * 100))
        coords = normalise_rows(coords)  # scale all rows to unit length
        p = self.kmeans(nclusters, coords)
        return p

    def mds_decomp(self, noise=False):

        matrix = self.get_dm(noise)

        dbc = double_centre(matrix)
        return eigen(dbc)

    def mds_cluster(
            self,
            nclusters,
            decomp,
            verbosity=0,
    ):

        l = np.diag(np.sqrt(np.abs(decomp.vals[:nclusters])))
        e = decomp.vecs[:, :nclusters]
        cve = decomp.cve
        coords = e.dot(l)
        if verbosity > 0:
            print('{0} dimensions explain {1:.2f}% of '
                  'the variance'.format(coords.shape[1], cve * 100))
        p = self.kmeans(nclusters, coords)
        return p

    @staticmethod
    def kmeans(nclusters, coords):
        est = KMeans(n_clusters=nclusters, n_init=50, max_iter=500)
        est.fit(coords)
        return Partition(est.labels_)

    @staticmethod
    def gmm(nclusters, coords, n_init=50, n_iter=500):
        est = GMM(n_components=nclusters, n_init=n_init, n_iter=n_iter)
        est.fit(coords)
        return Partition(est.predict(coords))

    # def plot_dendrogram(self, compound_key):
    #     """ Extracts data from clustering to plot dendrogram """
    #
    #     partition = self.partitions[compound_key]
    #     (linkmat, names, threshold) = self.plotting_info[compound_key]
    #     fig = plt.figure(figsize=(11.7, 8.3))
    #     dendrogram(
    #         linkmat,
    #         color_threshold=threshold,
    #         leaf_font_size=8,
    #         leaf_rotation=90,
    #         leaf_label_func=lambda leaf: names[leaf] + '_' + str(partition[leaf]),
    #         count_sort=True,
    #     )
    #     plt.suptitle('Dendrogram', fontsize=16)
    #     plt.title('Distance metric: {0}    Linkage method: {1}    Number of classes: {2}'.format(compound_key[0],
    #                                                                                              compound_key[1],
    #                                                                                              compound_key[2]),
    #               fontsize=12)
    #     plt.axhline(threshold, color='grey', ls='dashed')
    #     plt.xlabel('Gene')
    #     plt.ylabel('Distance')
    #     return fig
