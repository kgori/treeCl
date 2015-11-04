#!/usr/bin/env python
from __future__ import print_function

# standard library
import itertools

# third party
import matplotlib.pyplot as plt
from matplotlib import cm as CM
import numpy as np

# treeCl
from .collection import Collection
from .clustering import Spectral, MultidimensionalScaling
from .distance_matrix import DistanceMatrix
from .partition import Partition
from .errors import optioncheck
from .utils import flatten_list

import logging
logger = logging.getLogger(__name__)


class Plotter(object):
    def __init__(self, collection=None, records=None, dm=None,
                 metric='geo', **kwargs):

        """ Initialisation:
        A fully-formed Collection object can be given, or a list of records.
        If a list of records is provided a Collection object is built from it,
        with **kwargs being passed to the constructor.
        Additionally a distance matrix can be passed, or else one is constructed
        using the metric specified (default=geodesic)
        """

        if records:
            self.collection = Collection(records, **kwargs)
        else:
            self.collection = collection

        if dm is not None:
            self.dm = dm
        # else:
        #     self.dm = self.calc_dm(metric)

        self._warnings()

    def __len__(self):
        if self.collection:
            return len(self.collection.records)
        return 0

    def _warnings(self):
        if self.collection:
            if any([(r.tree is None) for r in self.collection.records]):
                logger.warn('No trees have been calculated for these records')

    # def calc_dm(self, method='geo'):
    #
    #     return (self.collection.distance_matrix(method)
    #             if self.collection
    #             else None)

    def get_decomp(self, method='MDS', **kwargs):
        optioncheck(method, ['MDS', 'spectral'])
        cl = Clustering(self.dm)
        if method == 'MDS':
            return cl.mds_decomp()
        if method == 'spectral':
            return cl.spectral_decomp(**kwargs)

    def get_coords(self, method, dimensions, normalise=False, **kwargs):
        decomp = self.get_decomp(method, **kwargs)
        if method in ['mds', 'MDS']:
            L = np.diag(np.sqrt(decomp.vals[:dimensions]))
            E = decomp.vecs[:, :dimensions]
            coords = E.dot(L)
            return coords

        else:
            return self.decomp_to_coords(decomp, dimensions, normalise=True)

    def decomp_to_coords(self, decomp, dimensions, normalise=False):
        optioncheck(dimensions, [2, 3])

        coords = decomp.coords_by_dimension(dimensions)[0]
        return coords.normalise_rows() if normalise else coords

    def sphere(self, ax):
        (u, v) = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        ax.plot_wireframe(x, y, z, color='grey', linewidth=0.2)
        return ax

    def heatmap(self, partition=None, cmap=CM.Blues):
        """ Plots a visual representation of a distance matrix """

        if isinstance(self.dm, DistanceMatrix):
            length = self.dm.values.shape[0]
        else:
            length = self.dm.shape[0]
        datamax = float(np.abs(self.dm).max())
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ticks_at = [0, 0.5 * datamax, datamax]
        if partition:
            sorting = flatten_list(partition.get_membership())
            self.dm = self.dm.reorder(sorting)
        cax = ax.imshow(
            self.dm.values,
            interpolation='nearest',
            origin='lower',
            extent=[0., length, 0., length],
            vmin=0,
            vmax=datamax,
            cmap=cmap,
        )
        cbar = fig.colorbar(cax, ticks=ticks_at, format='%1.2g')
        cbar.set_label('Distance')
        return fig

    def embedding(self, method='MDS', dimensions=3, partition=None,
                  add_sphere=False,
                  xlab='PCo1', ylab='PCo2', zlab='PCo3',
                  title='Trees embedded in dimension-reduced space',
                  outfile=False, **kwargs
    ):

        """ Gets coordinates, then calls embedding_plotter to do the plot """

        coords = self.get_coords(method, dimensions, normalise=False, **kwargs)
        return self.embedding_plotter(coords, dimensions, partition, add_sphere,
                                      xlab, ylab, zlab, title, outfile)

    def embedding_plotter(
            self, coordinates, dimensions, partition=None, add_sphere=False,
            xlab='PCo1', ylab='PCo2', zlab='PCo3',
            title='Trees embedded in dimension-reduced space',
            outfile=False,
    ):
        """ Points are coloured according to cluster membership specified
        by Partition object (or all black if no Partition specified) """

        optioncheck(dimensions, [2, 3])
        partition = (partition or
                     Partition(tuple([0] * len(coordinates))))

        colours = zip(*zip(range(len(partition)), itertools.cycle('bgrcmyk')))[1]
        print(colours)
        colour_mapping = np.array([colours[i - 1]
                                   for i in partition.partition_vector])
        fig = plt.figure()

        if dimensions == 3:
            ax = fig.add_subplot(111, projection='3d',
                                 xlabel=xlab, ylabel=ylab, zlabel=zlab, title=title)
            if add_sphere:
                ax = self.sphere(ax)

        else:
            ax = fig.add_subplot(111,
                                 xlabel=xlab, ylabel=ylab, title=title)

        ax.scatter(*coordinates.T, color=colour_mapping)
        # ax.set_aspect(1)


        if outfile:
            fig.savefig('{0}.pdf'.format(outfile))

        return fig


if __name__ == '__main__':
    # TESTS

    # from collection import Collection
    from lib.remote.utils import fileIO

    path_to_file = fileIO.path_to(__file__)
    test_data = fileIO.join_path(path_to_file, 'aa_alignments')

    print('Loading data...', end='')
    c = Collection(
        input_dir=test_data,
        file_format='phylip',
        datatype='protein',
        compression='gz',
    )
    print(' success')
    print('Building trees... ', end='')
    c.calc_NJ_trees()
    print('success')

    print('Building distance matrix...', end='')
    dm = c.distance_matrix('geo')
    print('success')

    print('Can build empty Plotter object... ', end='')
    empty_plotter = Plotter()
    print('true')

    print('Can build Plotter from Collection object... ', end='')
    plotter_from_collection = Plotter(c, metric='euc')
    print('yep')

    print('Can build Plotter from Collection + DistanceMatrix...', end='')
    plotter_with_dm = Plotter(c, dm=dm)
    print('yes')

    print('Can build Plotter from a list of TrClSeq objects...', end='')
    plotter_from_records = Plotter(records=c.records)
    print('yes')

    print('Can build Plotter from DistanceMatrix only...', end='')
    plotter_just_dm = Plotter(dm=dm)
    print('yes')

    print('Testing plotting')
    p = Partition(tuple([1] * 15 + [2] * 15 + [3] * 15 + [4] * 15))
    p_rand = Partition(tuple([1, 3, 1, 4, 2, 3, 3, 3, 2, 2, 1, 3, 3, 4, 1, 4, 1,
                              1, 2, 4, 1, 2, 2, 2, 2, 2, 3, 4, 2, 2, 1, 4, 3, 1, 4, 4, 3, 1, 3, 1, 3,
                              2, 4, 4, 1, 4, 1, 2, 3, 4, 2, 4, 3, 2, 1, 3, 4, 4, 1, 3]))
    fig1 = plotter_from_collection.embedding('MDS', 2, p)  # 2d MDS embedding
    fig2 = plotter_from_collection.embedding('MDS', 3, p)  # 3d MDS embedding
    fig3 = plotter_from_collection.embedding('spectral', 2, p_rand)  # 2d spectral
    fig4 = plotter_from_collection.embedding('spectral', 3, p_rand)  # 3d spectral
    fig5 = plotter_just_dm.heatmap(p)  # distance matrix as
    fig6 = plotter_just_dm.heatmap(p_rand)
    plt.show()
