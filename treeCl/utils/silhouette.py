from __future__ import print_function
from __future__ import division
from builtins import object
import numpy as np
import pandas as pd
from ..partition import Partition


class Silhouette(object):
    def __init__(self, dm):
        self._pvec = None
        self.distances = dm
        self.groups = None
        self.neighbours = None
        self.scores = None

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

    def get_indices_for_group(self, group):
        return np.where(self.pvec == group)[0]

    def get_indices_for_groups(self, group1, group2):
        ix = np.where(self.pvec == group1)[0]
        jx = np.where(self.pvec == group2)[0]
        return self.__get_indices_for_groups_by_index(ix, jx)

    def get_mean_dissimilarities_for_group(self, group):
        outgroups = self.groups[self.groups != group]
        within_indices = self.get_indices_for_groups(group, group)
        within_distances = self.distances[within_indices].mean(axis=1)
        dissimilarities = []
        for outgroup in outgroups:
            between_indices = self.get_indices_for_groups(group, outgroup)
            between_distances = self.distances[between_indices]
            dissimilarities.append(between_distances.mean(axis=1))
        return within_distances, np.array(dissimilarities), outgroups

    def run(self):
        if len(self.groups) == 1:
            raise ValueError("Silhouette is not defined for singleton clusters")
        for ingroup in self.groups:
            ingroup_ix = self.get_indices_for_group(ingroup)
            within, between, outgroups = self.get_mean_dissimilarities_for_group(ingroup)
            between_min = between.min(axis=0)
            outgroup_ix, neighbours_ix = np.where(between == between_min)
            neighbours = np.zeros(neighbours_ix.shape)
            neighbours[neighbours_ix] = outgroups[outgroup_ix]
            self.neighbours[ingroup_ix] = neighbours
            self.scores[ingroup_ix] = self.__silhouette_calc(within, between_min)

    @property
    def pvec(self):
        return self._pvec

    @pvec.setter
    def pvec(self, partition):
        if isinstance(partition, Partition):
            self._pvec = np.array(partition.partition_vector)
        else:
            self._pvec = np.array(partition)
        self.groups = np.unique(self._pvec)
        self.neighbours = np.zeros(self._pvec.shape)
        self.scores = np.zeros(self._pvec.shape)

    def __call__(self, partition):
        self.pvec = partition
        self.run()
        return self.neighbours, self.scores


def add_silhouettes_to_dataframe(path_to_distances, path_to_table, **kwargs):
    table = pd.read_csv(path_to_table, **kwargs)
    dm = np.loadtxt(path_to_distances)


if __name__ == '__main__':

    dm = np.array(
        [[0., 0.352, 0.23, 0.713, 0.426, 0.653, 0.481, 0.554, 1.533, 1.549, 1.505, 1.46],
         [0.352, 0., 0.249, 0.772, 0.625, 0.909, 0.668, 0.725, 1.613, 1.623, 1.568, 1.523],
         [0.23, 0.249, 0., 0.811, 0.417, 0.751, 0.456, 0.52, 1.489, 1.501, 1.446, 1.396],
         [0.713, 0.772, 0.811, 0., 0.962, 0.894, 1.025, 1.068, 1.748, 1.782, 1.724, 1.72],
         [0.426, 0.625, 0.417, 0.962, 0., 0.644, 0.083, 0.216, 1.424, 1.439, 1.398, 1.339],
         [0.653, 0.909, 0.751, 0.894, 0.644, 0., 0.685, 0.659, 1.467, 1.502, 1.448, 1.416],
         [0.481, 0.668, 0.456, 1.025, 0.083, 0.685, 0., 0.203, 1.419, 1.432, 1.394, 1.331],
         [0.554, 0.725, 0.52, 1.068, 0.216, 0.659, 0.203, 0., 1.503, 1.53, 1.472, 1.416],
         [1.533, 1.613, 1.489, 1.748, 1.424, 1.467, 1.419, 1.503, 0., 0.288, 0.299, 0.262],
         [1.549, 1.623, 1.501, 1.782, 1.439, 1.502, 1.432, 1.53, 0.288, 0., 0.296, 0.185],
         [1.505, 1.568, 1.446, 1.724, 1.398, 1.448, 1.394, 1.472, 0.299, 0.296, 0., 0.197],
         [1.46, 1.523, 1.396, 1.72, 1.339, 1.416, 1.331, 1.416, 0.262, 0.185, 0.197, 0.]])

    plist = [Partition((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)),
             Partition((1, 2, 3, 4, 5, 6, 5, 7, 8, 9, 10, 11)),
             Partition((1, 2, 3, 4, 5, 6, 5, 7, 8, 9, 10, 9)),
             Partition((1, 2, 1, 3, 4, 5, 4, 6, 7, 8, 9, 8)),
             Partition((1, 2, 1, 3, 4, 5, 4, 4, 6, 7, 8, 7)),
             Partition((1, 2, 1, 3, 4, 5, 4, 4, 6, 7, 7, 7)),
             Partition((1, 2, 1, 3, 4, 5, 4, 4, 6, 6, 6, 6)),
             Partition((1, 1, 1, 2, 3, 4, 3, 3, 5, 5, 5, 5)),
             Partition((1, 1, 1, 2, 3, 3, 3, 3, 4, 4, 4, 4)),
             Partition((1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)),
             Partition((1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2)),
             Partition((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))]

    s = Silhouette(dm)
    skips = 0
    for p in plist:

        try:
            neighbours, scores = s(p)
            print("{} clusters: avg score = {}".format(len(p), scores.mean()))
        except ValueError:
            print("{} clusters: skipping".format(len(p)))
            skips += 1
    print ("{} tests, {} skipped".format(len(plist), skips))

    import treeCl
    cl = treeCl.Clustering(dm)
    skips = 0
    for p in plist:
        try:
            anosim = cl.anosim(p)
        except ValueError:
            skips += 1
            continue
        try:
            permanova = cl.permanova(p)
        except ValueError:
            skips += 1
            continue
        print ("{} clusters: anosim = {}; permanova = {}".format(len(p), anosim.p_value, permanova.p_value))
    print ("{} tests, {} skipped".format(2*len(plist), skips))

