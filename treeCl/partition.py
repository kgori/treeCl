from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
__author__ = 'kgori'

from collections import defaultdict
import numpy as np
import random
import scipy.stats

def entropies(partition_1, partition_2):
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

    if partition_1.num_elements() != partition_2.num_elements():
        print('Partition lists are not the same length')
        return 0
    else:
        total = partition_1.num_elements()

    m1 = partition_1.get_membership()
    m2 = partition_2.get_membership()
    l1 = len(m1)
    l2 = len(m2)
    entropy_1 = 0
    entropy_2 = 0
    mut_inf = 0
    for i in range(l1):
        prob1 = len(m1[i]) / total  # float division ensured by __future__ import
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

    return entropy_1, entropy_2, mut_inf


class Partition(object):
    """ Class to store clustering information """

    def __init__(self, partition_vector):
        self._partition_vector = None
        self.partition_vector = partition_vector
        self.membership = self.get_membership()

    def __str__(self):
        return str(self.partition_vector)

    def __repr__(self):
        return self.__class__.__name__ + '({0})'.format(str(self))

    def __len__(self):
        """
        This gives the number of groups in the partition, rather than
        the number of elements in the data
        """
        return max(self.partition_vector) + 1

    def __getitem__(self, index):
        return self.membership[index]

    @property
    def partition_vector(self):
        return self._partition_vector

    @partition_vector.setter
    def partition_vector(self, vec):
        self._partition_vector = self._restricted_growth_notation(vec)

    def is_minimal(self):
        """ The partition describes all members being in the same cluster
        :return: boolean
        """
        return self.num_groups() == 1

    def is_maximal(self):
        """ The partition describes every member being in its own group
        :return: boolean
        """
        return self.num_groups() == self.num_elements()

    @staticmethod
    def _restricted_growth_notation(l):
        """ The clustering returned by the hcluster module gives group
        membership without regard for numerical order This function preserves
        the group membership, but sorts the labelling into numerical order """

        list_length = len(l)

        d = defaultdict(list)
        for (i, element) in enumerate(l):
            d[element].append(i)

        l2 = [None] * list_length

        for (name, index_list) in enumerate(sorted(d.values(), key=min)):
            for index in index_list:
                l2[index] = name

        return tuple(l2)

    @classmethod
    def random(cls, alpha, size):
        """
        Generate a random start using expected proportions, alpha.
        These are used to parameterise a random draw from a Dirichlet
        distribution.
        An example, to split a dataset of 20 items into 3 groups of [10,
        6, 4] items:
         - alpha = [10, 6, 4],
         - alpha = [100, 60, 40],
         - alpha = [5, 3, 2],
        would all work. Variance is inversely related to sum(alpha)
        """
        props = np.concatenate([[0], (scipy.stats.dirichlet.rvs(alpha) * size).cumsum().round().astype(int)])
        indices = np.array(list(range(size)))
        random.shuffle(indices)
        x = []
        for i in range(len(props)-1):
            ix = indices[props[i]:props[i+1]]
            x.append(ix)
        return cls.from_membership(x)

    def num_elements(self):
        return len(self.partition_vector)

    def num_groups(self):
        return max(self.partition_vector) + 1

    @classmethod
    def from_membership(cls, membership):
        tmp = {}
        for group, members in enumerate(membership):
            for member in members:
                tmp[member] = group
        return cls([tmp[k] for k in sorted(tmp)])

    def get_membership(self):
        """
        Alternative representation of group membership -
        creates a list with one tuple per group; each tuple contains
        the indices of its members

        Example:
        partition  = (0,0,0,1,0,1,2,2)
        membership = [(0,1,2,4), (3,5), (6,7)]

        :return: list of tuples giving group memberships by index
        """
        result = defaultdict(list)
        for (position, value) in enumerate(self.partition_vector):
            result[value].append(position)
        return sorted([tuple(x) for x in result.values()])

    @classmethod
    def read(cls, filename):
        with open(filename) as reader:
            s = reader.read()
            t = tuple(s.rstrip().split(','))
        return cls(t)

    def normalised_mutual_information(self, other):
        (entropy_1, entropy_2, mut_inf) = entropies(self, other)

        return 2 * mut_inf / (entropy_1 + entropy_2)

    def variation_of_information(self, other):
        """ calculates Variation of Information Metric between two clusterings
        of the same data - SEE Meila, M. (2007). Comparing clusterings: an
        information based distance. Journal of Multivariate Analysis, 98(5),
        873-895. doi:10.1016/j.jmva.2006.11.013"""
        (entropy_1, entropy_2, mut_inf) = entropies(self, other)

        return entropy_1 + entropy_2 - 2 * mut_inf

    def write(self, filename):
        with open(filename, 'w') as writer:
            writer.write(','.join(str(x) for x in self.partition_vector) + '\n')
