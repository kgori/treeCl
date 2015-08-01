#!/usr/bin/env python
from __future__ import print_function

# standard library

# third party
import numpy as np
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
import fastcluster
import skbio

try:
    from Bio.Cluster import kmedoids
    Biopython_Unavailable = False
except ImportError:
    print("Biopython unavailable - kmedoids clustering disabled")
    Biopython_Unavailable = True

from sklearn.cluster import AffinityPropagation, DBSCAN, KMeans
from sklearn.mixture import GMM
from sklearn.manifold import spectral_embedding

# treeCl
from .distance_matrix import DistanceMatrix, rbf, binsearch_mask, kmask, kscale, affinity, laplace, eigen, double_centre, \
    normalise_rows
from .partition import Partition
from .utils import enum
from .errors import OptionError, isnumbercheck, rangecheck

options = enum(
    "PRUNING_NONE",
    "PRUNING_ESTIMATE",
    "PRUNING_MANUAL",
    "LOCAL_SCALE_MEDIAN",
    "LOCAL_SCALE_ESTIMATE",
    "LOCAL_SCALE_MANUAL")

methods = enum(
    "KMEANS",
    "GMM")

linkage = enum(
    "SINGLE",
    "COMPLETE",
    "AVERAGE",
    "WARD",
    "WEIGHTED",
    "CENTROID",
    "MEDIAN")


class ClusteringManager(object):
    """
    Clustering manager base class
    """
    def __init__(self, dm):
        if isinstance(dm, np.ndarray):
            dm = DistanceMatrix.from_array(dm)

        if not isinstance(dm, DistanceMatrix):
            raise ValueError('Distance matrix should be a numpy array or treeCl.DistanceMatrix')
        self.dm = dm

    def __str__(self):
        return str(self.dm)

    def get_dm(self, noise):
        return self.dm.add_noise().values if noise else self.dm.values


class EMMixin(object):
    """
    Provide methods to do kmeans and GMM estimation
    """
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


def _check_val(opt, min_, max_):
    isnumbercheck(opt)
    rangecheck(opt, min_, max_)

class Spectral(ClusteringManager, EMMixin):
    """
    Manager for spectral clustering 
    """
    def __init__(self, dm, 
                 pruning_option=options.PRUNING_NONE, 
                 scale_option=options.LOCAL_SCALE_MEDIAN, 
                 manual_pruning=None, 
                 manual_scale=None,
                 verbosity=0):
        super(Spectral, self).__init__(dm)
        try:
            options.reverse[pruning_option]
        except KeyError:
            raise OptionError(pruning_option, options.reverse.values())
        try:
            options.reverse[scale_option]
        except KeyError:
            raise OptionError(scale_option, options.reverse.values())

        if pruning_option == options.PRUNING_MANUAL:
            _check_val(manual_pruning, 2, self.dm.df.shape[0])

        if scale_option == options.LOCAL_SCALE_MANUAL:
            _check_val(manual_scale, 2, self.dm.df.shape[0])

        self._pruning_option = pruning_option
        self._scale_option = scale_option
        self._manual_pruning = manual_pruning
        self._manual_scale = manual_scale
        self._verbosity = verbosity
        self._affinity = self.decompose()

    def __str__(self):
        return ('Spectral Clustering with local scaling:\n'
                'Pruning option: {}\n'
                'Scaling option: {}'
                .format(options.reverse[self._pruning_option], options.reverse[self._scale_option]))

    def decompose(self,
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

        # get local scale estimate
        est_scale = None

        # ADJUST MASK
        if self._pruning_option == options.PRUNING_NONE:  
            # Set kp to max value
            kp = len(matrix) - 1
            mask = np.ones(matrix.shape, dtype=bool)
        elif self._pruning_option == options.PRUNING_MANUAL:
            # Manually set value of kp
            kp = self._manual_pruning
            mask = kmask(matrix, self._manual_pruning, logic=logic)
        elif self._pruning_option == options.PRUNING_ESTIMATE:
            # Must estimate value of kp
            kp, mask, est_scale = binsearch_mask(matrix, logic=logic) 
        else:
            raise ValueError("Unexpected error: 'kp' not set")

        # ADJUST SCALE
        if self._scale_option == options.LOCAL_SCALE_MEDIAN:
            dist = np.median(matrix, axis=1)
            scale = np.outer(dist, dist)
        elif self._scale_option == options.LOCAL_SCALE_MANUAL:
            scale = kscale(matrix, self._manual_scale)
        elif self._scale_option == options.LOCAL_SCALE_ESTIMATE:
            if est_scale is None:
                _, _, scale = binsearch_mask(matrix, logic=logic) 
            else:
                # Nothing to be done - est_scale was set during the PRUNING_ESTIMATE
                scale = est_scale
        else:
            raise ValueError("Unexpected error: 'scale' not set")

        # ZeroDivisionError safety check
        if not (scale > 1e-5).all():
            if verbosity > 0:
                print('Rescaling to avoid zero-div error')
            _, _, scale = binsearch_mask(matrix, logic=logic)
            assert (scale > 1e-5).all()

        aff = affinity(matrix, mask, scale)
        return aff

    def cluster(self, n, method=methods.KMEANS):
        """
        Cluster the embedded coordinates

        Parameters
        ----------
        n:      int
                The number of clusters to return
        method: enum value, one of (methods.KMEANS | methods.GMM)
                The clustering method to use

        Returns
        -------
        Partition: Partition object describing the data partition
        """
        if n == 1:
            return Partition([1] * len(self.get_dm(False)))

        coords = spectral_embedding(self._affinity, n)
        self._coords = normalise_rows(coords)  # scale all rows to unit length
        if method == methods.KMEANS:
            p = self.kmeans(n, self._coords)
        elif method == methods.GMM:
            p = self.gmm(n, self._coords)
        else:
            raise OptionError(method, list(methods.reverse.values()))
        if self._verbosity > 0:
            print('Using clustering method: {}'.format(methods.reverse[method]))
        return p

    @property
    def affinity(self):
        return self._affinity
    


class Hierarchical(ClusteringManager):
    """ Apply clustering methods to distance matrix

    = Hierarchical clustering - single-linkage - complete-linkage - average-
    linkage (UPGMA) - Ward's method

    = k-medoids

    = Multidimensional Scaling (Principal Coordinate Analysis) + k-means

    = Spectral Clustering + k-means - NJW method - Shi-Malik method - Zelnik-
    Manor and Perona Local Scaling - Local Scaling with eigenvector rotation as
    stopping criterion

    """

    def __str__(self):
        return 'Hierarchical Clustering'

    def cluster(self, nclusters, linkage_method=linkage.WARD, **kwargs):
        """
        Do hierarchical clustering on a distance matrix using one of the methods:
            methods.SINGLE   = single-linkage clustering
            methods.COMPLETE = complete-linkage clustering
            methods.AVERAGE  = average-linkage clustering
            methods.WARD     = Ward's minimum variance method
        """

        if linkage_method == linkage.SINGLE:
            return self._hclust(nclusters, 'single', **kwargs)
        elif linkage_method == linkage.COMPLETE:
            return self._hclust(nclusters, 'complete', **kwargs)
        elif linkage_method == linkage.AVERAGE:
            return self._hclust(nclusters, 'average', **kwargs)
        elif linkage_method == linkage.WARD:
            return self._hclust(nclusters, 'ward', **kwargs)
        elif linkage_method == linkage.WEIGHTED:
            return self._hclust(nclusters, 'weighted', **kwargs)
        elif linkage_method == linkage.CENTROID:
            return self._hclust(nclusters, 'centroid', **kwargs)
        elif linkage_method == linkage.MEDIAN:
            return self._hclust(nclusters, 'median', **kwargs)
        else:
            raise ValueError('Unknown linkage_method: {}'.format(linkage_method))

    def _hclust(self, nclusters, method, noise=False):
        """
        :param nclusters: Number of clusters to return
        :param linkage_method: single, complete, average, ward, weighted, centroid or median
                               (http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html)
        :param noise: Add Gaussian noise to the distance matrix prior to clustering (bool, default=False)
        :return: Partition object describing clustering
        """
        matrix = self.get_dm(noise)

        linkmat = fastcluster.linkage(squareform(matrix), method)
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


class Clustering(ClusteringManager, EMMixin):
    """ Apply clustering methods to distance matrix

    = Hierarchical clustering - single-linkage - complete-linkage - average-
    linkage (UPGMA) - Ward's method

    = k-medoids

    = Multidimensional Scaling (Principal Coordinate Analysis) + k-means

    = Spectral Clustering + k-means - NJW method - Shi-Malik method - Zelnik-
    Manor and Perona Local Scaling - Local Scaling with eigenvector rotation as
    stopping criterion

    """

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


class Assessment(ClusteringManager):

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

    def silhouette(self, partition):
        pvec   = np.array(partition.partition_vector)
        groups = np.unique(pvec)
        nbrs   = np.zeros(pvec.shape)
        scores = np.zeros(pvec.shape)

        if len(groups) == 1:
            raise ValueError("Silhouette is not defined for singleton clusters")
        for ingroup in groups:
            ingroup_ix = np.where(pvec == ingroup)[0]
            within, between, outgroups = self.__get_mean_dissimilarities_for_group(pvec, ingroup, groups)
            between_min = between.min(axis=0)
            outgroup_ix, neighbours_ix = np.where(between == between_min)
            neighbours = np.zeros(neighbours_ix.shape)
            neighbours[neighbours_ix] = outgroups[outgroup_ix]
            nbrs[ingroup_ix] = neighbours
            scores[ingroup_ix] = self.__silhouette_calc(within, between_min)

        return nbrs, scores

    @staticmethod
    def __get_indices_for_groups_by_index(ix, jx):
        if len(ix) == len(jx) == 1 and ix == jx:
            return [list(ix)], [list(jx)]
        row_indices = [[i for j in jx if i != j] for i in ix]
        column_indices = [[j for j in jx if j != i] for i in ix]
        return row_indices, column_indices

    @staticmethod
    def __silhouette_calc(ingroup, outgroup):
        if len(ingroup) == 1:
            return 0
        max_ = np.array([ingroup, outgroup]).max(axis=0)
        return (outgroup - ingroup) / max_

    def __get_indices_for_groups(self, pvec, group1, group2):
        ix = np.where(pvec == group1)[0]
        jx = np.where(pvec == group2)[0]
        return self.__get_indices_for_groups_by_index(ix, jx)

    def __get_mean_dissimilarities_for_group(self, pvec, group, groups):
        outgroups = groups[groups != group]
        within_indices = self.__get_indices_for_groups(pvec, group, group)
        within_distances = self.dm.values[within_indices].mean(axis=1)
        dissimilarities = []
        for outgroup in outgroups:
            between_indices = self.__get_indices_for_groups(pvec, group, outgroup)
            between_distances = self.dm.values[between_indices]
            dissimilarities.append(between_distances.mean(axis=1))
        return within_distances, np.array(dissimilarities), outgroups



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
