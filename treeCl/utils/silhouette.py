import numpy as np
import pandas as pd


class Silhouette(object):
    def __init__(self, dm, p):
        self._partition = None
        self.partition = p
        self.distances = dm
        self.groups = None
        self.neighbours = None
        self.scores = None

    def get_indices_for_group(self, group):
        return np.where(self.partition == group)[0]

    def get_indices_for_groups(self, group1, group2):
        ix = np.where(self.partition == group1)[0]
        jx = np.where(self.partition == group2)[0]
        return self.__get_indices_for_groups_by_index(ix, jx)

    def __get_indices_for_groups_by_index(self, ix, jx):
        if len(ix) == len(jx) == 1 and ix == jx:
            return [list(ix)], [list(jx)]
        row_indices = [[i for j in jx if i != j] for i in ix]
        column_indices = [[j for j in jx if j != i] for i in ix]
        return row_indices, column_indices

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

    def __silhouette_calc(self, ingroup, outgroup):
        if len(ingroup) == 1:
            return 0
        max_ = np.array([ingroup, outgroup]).max(axis=0)
        return (outgroup - ingroup) / max_

    def run(self):
        for ingroup in self.groups:
            ingroup_ix = self.get_indices_for_group(ingroup)
            within, between, outgroups = self.get_mean_dissimilarities_for_group(ingroup)
            between_min = between.min(axis=0)
            outgroup_ix, neighbours_ix = np.where(between == between_min)
            neighbours = np.zeros(neighbours_ix.shape)
            neighbours[neighbours_ix] = outgroups[outgroup_ix]
            self.neighbours[ingroup_ix] = neighbours
            self.scores[ingroup_ix] = self.__silhouette_calc(within,
                                                             between_min)

    @property
    def partition(self):
        return self._partition

    @partition.setter
    def partition(self, partition):
        if isinstance(partition, treeCl.Partition):
            self.partition = np.array(partition.partition_vector)
        else:
            self.partition = np.array(partition)
        self.groups = np.unique(self._partition)
        self.neighbours = np.zeros(self._partition.shape)
        self.scores = np.zeros(self._partition.shape)

    def silhouette(self, partition):
        self.partition = partition
        self.run()
        return (self.neighbours, self.scores)


def add_silhouettes_to_dataframe(path_to_distances, path_to_table, **kwargs):
    table = pd.read_csv(path_to_table, **kwargs)
    dm = np.loadtxt(path_to_distances)


s = Silhouette(dm, p)
s.run()
