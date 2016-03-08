from __future__ import division, print_function
from builtins import zip
from builtins import str
from builtins import object

import numpy as np

from ..partition import Partition


__author__ = 'kgori'


class Kelley(object):
    def __init__(self, dm):
        self._dm = None
        self.dm = dm
        self.n_obs = dm.shape[0]

    @property
    def dm(self):
        return self._dm

    @dm.setter
    def dm(self, dm):
        self._dm = dm

    def _average_spread(self, partition):
        # print('_as() partition: {}'.format(partition))
        indices = self.partition_to_indices(partition)
        # print('_as() indices: {}'.format(indices))
        avg_spread = np.mean([self._spread(*ix) for ix in indices])
        # print('_as() result: {}'.format(avg_spread))
        return avg_spread

    def _normalised_average_spread(self, *partitions):
        partitions = [pt for pt in partitions if 1 < len(pt) < self.n_obs]
        # print('_nas()', partitions)
        avg_spread_list = np.array([self._average_spread(partition) for partition in partitions])
        # print('_nas() avg_spread_list: {}'.format(avg_spread_list))
        if len(partitions) == 1:
            return avg_spread_list
        max_ = avg_spread_list.max()
        min_ = avg_spread_list.min()
        range_ = max_ - min_
        N = self.dm.shape[0]
        avg_spread_norm = ((N - 2) / range_) * (avg_spread_list - min_) + 1
        # print('_nas() result {}'.format(avg_spread_norm))
        return avg_spread_norm

    def partition_to_indices(self, partition):
        # print('_pti', partition)
        ix = [tup for tup in partition.get_membership() if len(tup) > 1]
        return ix

    def _spread(self, *ix):
        # print('_spread() ix: {}'.format(ix))
        N = len(ix)
        if N == 1:
            result = 0.
        else:
            ix = np.array(ix)[np.newaxis]
            cluster_m = self.dm[ix, ix.T]
            # print(cluster_m)
            upper_tri = cluster_m[np.triu_indices(N, 1)]
            result = upper_tri.sum() / (0.5 * self.n_obs * (self.n_obs - 1))
        # print('_spread() result: {}'.format(result))
        return result

    def penalty_values(self, *partitions):
        partitions = [p for p in partitions if 1 < len(p) < self.n_obs]
        # print('pv()', partitions)
        assert isinstance(partitions, (list, tuple))
        avg_spread_norm = self._normalised_average_spread(*partitions)
        nC = (len(p) for p in partitions)
        return dict(zip(nC, avg_spread_norm + np.array([len(p) for p in partitions])))

    def pen_val_debug(self, *partitions):
        partitions = [p for p in partitions if 1 < len(p) < self.n_obs]
        nc = np.array([len(p) for p in partitions])
        print('clusters per partition: {}'.format([str(x) for x in nc]))
        avg_spread = np.array([self._average_spread(partition) for partition in partitions])
        print('average_spread_list: {}'.format(list(avg_spread)))
        avg_spread_norm = self._normalised_average_spread(*partitions)
        print('normalised: {}'.format(avg_spread_norm))
        return dict(zip(nc, avg_spread_norm + nc))


if __name__ == '__main__':

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

    k = Kelley(dm)
    for p in plist:
        print(k.partition_to_indices(p))
    print(k.penalty_values(*plist))
    # print(k._spread(8,9,10,11))
    # print(k.pen_val_debug(*plist[0:11]))
    # print(k.pen_val_debug(*plist[1:12]))
    # print(k.pen_val_debug(*plist[1:11]))
