#!/usr/bin/env python

import numpy as np
from distance_matrix import DistanceMatrix
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
try:
    from Bio.Cluster import kmedoids
    Biopython_Unavailable = False
except ImportError:
    print "Biopython unavailable - kmedoids clustering disabled"
    Biopython_Unavailable = True
import matplotlib.pyplot as plt
try:
    from sklearn.cluster import KMeans
except ImportError:
    print "sklearn unavailable: KMeans unavailable"
from collections import defaultdict
#import evrot ## evrot not currently in use
from copy import deepcopy


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

    def __init__(self, distance_matrix):

        self.distance_matrix = distance_matrix

    def __str__(self):
        return str(self.distance_matrix)

    def kmedoids(self, nclusters, noise=False):

        if Biopython_Unavailable:
            print 'kmedoids not available without Biopython'
            return
        
        if noise:
            matrix = self.distance_matrix.add_noise()
        else:
            matrix = self.distance_matrix

        p = [kmedoids(matrix, nclusters=nclusters, npass=100) for _ in
             range(100)]
        p.sort(key=lambda x: x[1])
        return Partition(p[0][0])

    def hierarchical(
        self,
        nclusters,
        linkage_method,
        noise=False,
        ):

        if noise:
            matrix = self.distance_matrix.add_noise()
        else:
            matrix = self.distance_matrix

        linkmat = linkage(matrix, linkage_method)
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
        T = fcluster(linkmat, threshold, criterion='distance')
        return Partition(T)

    def spectral_decomp(
        self,
        prune,
        local_scale=7,
        sigma=2,
        noise=False,
        ):
        """ Use prune to remove links between distant points:
        
        prune='estimate' searches for the smallest value that retains a fully
        connected graph
        
        """

        if noise:
            matrix = self.distance_matrix.add_noise()
        else:
            matrix = self.distance_matrix
        kp = None

        if prune == 'estimate':       # use binary search
            (kp, mask) = matrix.binsearch_mask()
            print 'Pruning distances greater than {0} nearest-neighbour'.format(kp)
        
        elif prune == -1:             # no pruning
            mask = None
            print 'No pruning applied'
        
        else:                         # prune to chosen value
            kp = prune
            mask = matrix.kmask(k=kp)
            print 'Pruning connections beyond {0} nearest-neighbour'.format(kp)

        if local_scale == 'estimate': # use same value as binary mask search 
            if kp is None: 
                kp = matrix.binsearch_mask()[0]
            scale = matrix.kdists(k=kp)
            print 'Local scale based on {0} nearest-neighbour'.format(kp)
        
        elif local_scale == 'median': # use median vector
            scale = np.median(matrix, axis=1)
        
        elif local_scale == -1:       # use maximum vector
            scale = matrix.kdists(k=matrix.shape[0])
            print 'Local scale based on {0} nearest-neighbour'.format(matrix.shape[0])
        
        else:                         # use chosen value
            scale = matrix.kdists(k=local_scale)
            print 'Local scale based on {0} nearest-neighbour'.format(local_scale)

        if not (scale > 0).all():
            (_, scale) = matrix.binsearch_dists()
        assert (scale > 0).all() 
        aff = matrix.affinity(mask, scale)
        laplace = matrix.laplace(aff)
        return laplace.eigen()  # vectors are in columns

    def spectral_cluster(self, nclusters, decomp):

        (coords, cve) = decomp.coords_by_dimension(nclusters)
        print '{0} dimensions explain {1:.2f}% of the variance'.format(nclusters,
                cve * 100)
        coords = coords.normalise_rows()  # scale all rows to unit length
        P = self.kmeans(nclusters, coords)
        return P

    def MDS_decomp(self, noise=False):

        if noise:
            matrix = self.distance_matrix.add_noise()
        else:
            matrix = self.distance_matrix

        dbc = matrix.double_centre()
        return dbc.eigen()

    def MDS_cluster(
        self,
        nclusters,
        decomp,
        cutoff=None,
        ):

        (coords, cve) = \
            (decomp.coords_by_cutoff(cutoff) if cutoff else decomp.coords_by_dimension(nclusters))
        print '{0} dimensions explain {1:.2f}% of the variance'.format(coords.shape[1],
                cve * 100)
        P = self.kmeans(nclusters, coords)
        return P

    def kmeans(self, nclusters, coords=None, noise=False):
        coords = coords if coords is not None else self.distance_matrix
        if noise:
            coords.add_noise()
        est = KMeans(n_clusters=nclusters, n_init=50, max_iter=500)
        est.fit(coords)
        return Partition(est.labels_)

    # NOT FIXED YET (26/02/13)

    def plot_dendrogram(self, compound_key):
        """ Extracts data from clustering to plot dendrogram """

        partition = self.partitions[compound_key]
        (linkmat, names, threshold) = self.plotting_info[compound_key]
        fig = plt.figure(figsize=(11.7, 8.3))
        dendrogram(
            linkmat,
            color_threshold=threshold,
            leaf_font_size=8,
            leaf_rotation=90,
            leaf_label_func=lambda leaf: names[leaf] + '_' \
                + str(partition[leaf]),
            count_sort=True,
            )
        plt.suptitle('Dendrogram', fontsize=16)
        plt.title('Distance metric: {0}    Linkage method: {1}    Number of classes: {2}'.format(compound_key[0],
                  compound_key[1], compound_key[2]), fontsize=12)
        plt.axhline(threshold, color='grey', ls='dashed')
        plt.xlabel('Gene')
        plt.ylabel('Distance')
        return fig

    # NOT FIXED YET (26/02/13)

    def spectral_rotate(
        self,
        prune=True,
        KMeans=True,
        recalculate=False,
        max_groups=None,
        min_groups=2,
        verbose=True,
        ):

        if dm.metric == 'rf':
            noise = True
        else:
            noise = False
        if recalculate or not 'spectral_decomp' in self.cache:
            laplacian = self.spectral(dm, prune=prune, add_noise=noise)

            (eigvals, eigvecs, cve) = self.get_eigen(laplacian,
                    standardize=False)
            self.cache['spectral_decomp'] = (eigvals, eigvecs, cve)
            self.cache['laplacian'] = laplacian
        else:

            (eigvals, eigvecs, cve) = self.cache['spectral_decomp']

        # ######################
        # CLUSTER_ROTATE STUFF HERE

        M = dm.matrix
        if not max_groups:
            max_groups = int(np.sqrt(M.shape[0]) + np.power(M.shape[0], 1.0
                             / 3))
        (nclusters, clustering, quality_scores, rotated_vectors) = \
            self.cluster_rotate(eigvecs, max_groups=max_groups,
                                min_groups=min_groups)

        translate_clustering = [None] * len(M)
        no_of_empty_clusters = 0
        for (group_number, group_membership) in enumerate(clustering):
            if len(group_membership) == 0:
                no_of_empty_clusters += 1
            for index in group_membership:
                translate_clustering[index - 1] = group_number
        clustering = Partition(translate_clustering)
        if no_of_empty_clusters > 0:
            print 'Subtracting {0} empty {1}'.format(no_of_empty_clusters,
                    ('cluster' if no_of_empty_clusters == 1 else 'clusters'))
            nclusters -= no_of_empty_clusters

        # ######################

        if verbose:
            print 'Discovered {0} clusters'.format(nclusters)
            print 'Quality scores: {0}'.format(quality_scores)
            if KMeans:
                print 'Pre-KMeans clustering: {0}'.format(clustering)
        if KMeans:
            T = self.kmeans(nclusters, rotated_vectors)
        else:
            T = clustering
        return (T, nclusters, quality_scores)

        # ######################

    # NOT FIXED YET (26.02.13)

    def cluster_rotate(
        self,
        eigenvectors,
        max_groups,
        min_groups=2,
        ):

        groups = range(min_groups, max_groups + 1)
        vector_length = eigenvectors.shape[0]
        current_vector = eigenvectors[:, :groups[0]]
        n = max_groups - min_groups + 1

        quality_scores = [None] * n
        clusters = [None] * n
        rotated_vectors = [None] * n

        for g in range(n):
            if g > 0:
                current_vector = np.concatenate((rotated_vectors[g - 1],
                        eigenvectors[:, groups[g] - 1:groups[g]]), axis=1)

            (clusters[g], quality_scores[g], rotated_vectors[g]) = \
                evrot.main(current_vector)

        # Find the highest index of quality scores where the
        # score is within 0.0025 of the maximum:
        # this is our chosen number of groups

        max_score = max(quality_scores)
        index = quality_scores.index(max_score)
        start = index + 1
        for (i, score) in enumerate(quality_scores[index + 1:], start=start):
            if abs(score - max_score) < 0.0025:
                index = i

        return (groups[index], clusters[index], quality_scores,
                rotated_vectors[index])


class Partition(object):

    """ Class to store clustering information """

    score = 0
    concats = []
    partition_vector = None

    def __init__(self, partition_vector):
        self.partition_vector = partition_vector

    def __str__(self):
        return str(self.partition_vector)

    def __repr__(self):
        return self.__class__.__name__ + str(self)

    def __len__(self):
        return len(self.partition_vector)

    @property
    def partition_vector(self):
        return self._partition_vector

    @partition_vector.setter
    def partition_vector(self, vec):
        self._partition_vector = self.order(vec)

    @classmethod
    def order(self, l):
        """ The clustering returned by the hcluster module gives group
        membership without regard for numerical order This function preserves
        the group membership, but sorts the labelling into numerical order """

        list_length = len(l)

        d = defaultdict(list)
        for (i, element) in enumerate(l):
            d[element].append(i)

        l2 = [None] * list_length

        for (name, index_list) in enumerate(sorted(d.values(), key=min),
                start=1):
            for index in index_list:
                l2[index] = name

        return tuple(l2)

    def entropies(self, partition_1, partition_2):
        """ parameters: partition_1 (list / array) - a partitioning of a dataset
        according to some clustering method. Cluster labels are arbitrary.
        partition_2 (list / array) - another partitioning of the same dataset.
        Labels don't need to match, nor do the number of clusters.
        
        subfunctions: get_membership( parameter partition (list / array) )
        returns a list of length equal to the number of clusters found in the
        partition. Each element is the set of members of the cluster. Ordering
        is arbitrary.
        
        variables used: t = total number of points in the dataset m1 = cluster
        memberships from partition_1 m2 = cluster memberships from partition_2
        l1 = length (i.e. number of clusters) of m1 l2 = length of m2 entropy_1
        = Shannon entropy of partition_1 entropy_2 = Shannon entropy of
        partition_2 mut_inf = mutual information of partitions prob1 =
        probability distribution of partition 1 - i.e. the probability that a
        randomly chosen datapoint belongs to each cluster (size of cluster /
        size of dataset) prob2 = as above, for partition 2 intersect = number of
        common elements in partition 1 [i] and partition 2 [j] """

        if len(partition_1) != len(partition_2):
            print 'Partition lists are not the same length'
            return 0
        else:
            total = float(len(partition_1))  # Ensure float division later

        m1 = self.get_membership(partition_1)
        m2 = self.get_membership(partition_2)
        l1 = len(m1)
        l2 = len(m2)
        entropy_1 = 0
        entropy_2 = 0
        mut_inf = 0
        for i in range(l1):
            prob1 = len(m1[i]) / total
            entropy_1 -= prob1 * np.log2(prob1)
            for j in range(l2):
                if i == 0:  # only calculate these once
                    prob2 = len(m2[j]) / total
                    entropy_2 -= prob2 * np.log2(prob2)
                intersect = len(set(m1[i]) & set(m2[j]))
                if intersect == 0:
                    continue  # because 0 * log(0) = 0 (lim x->0: xlog(x)->0)
                else:
                    mut_inf += intersect / total * np.log2(total * intersect
                            / (len(m1[i]) * len(m2[j])))

        return (entropy_1, entropy_2, mut_inf)

    def flatten(self, list_of_lists):
        """ This is faster than the one-liner version:-
        
        def(flatten): return list(itertools.chain(*list_of_lists))
        
        """

        res = []
        x = res.extend
        for sublist in list_of_lists:
            x(sublist)
        return res

    def get_membership(self, partition_vector=None, flatten=False):

        pvec = partition_vector or self.partition_vector
        result = defaultdict(list)
        for (position, value) in enumerate(pvec):
            result[value].append(position)
        result = [tuple(x) for x in sorted(result.values(), key=len,
                  reverse=True)]
        return (self.flatten(result) if flatten else result)

    def normalised_mutual_information(self, other):
        partition_1 = self.partition_vector
        partition_2 = other.partition_vector    

        (entropy_1, entropy_2, mut_inf) = self.entropies(partition_1,
                partition_2)

        return 2 * mut_inf / (entropy_1 + entropy_2)

    def variation_of_information(self, other):
        """ calculates Variation of Information Metric between two clusterings
        of the same data - SEE Meila, M. (2007). Comparing clusterings: an
        information based distance. Journal of Multivariate Analysis, 98(5),
        873-895. doi:10.1016/j.jmva.2006.11.013"""
        partition_1 = self.partition_vector
        partition_2 = other.partition_vector

        (entropy_1, entropy_2, mut_inf) = self.entropies(partition_1,
                partition_2)

        return entropy_1 + entropy_2 - 2 * mut_inf
