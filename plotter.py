#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import cm as CM
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from lib.remote.errors import optioncheck

class Plotter(object):

    def __init__(self):
        pass

    def decomp_to_coords(self, decomp, dimensions, normalise=False):
        optioncheck(dimensions, [2,3])

        coords = decomp.coords_by_dimension(dimensions)[0]
        return coords.normalise_rows() if normalise else coords

    def sphere(self, ax):
        (u, v) = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
        x = np.cos(u) * np.sin(v)
        y = np.sin(u) * np.sin(v)
        z = np.cos(v)
        ax.plot_wireframe(x, y, z, color='grey', linewidth=0.2)
        return ax

    def embedding(
            self, coordinates, dimensions, partition, add_sphere=False, 
            xlab='PCo1', ylab='PCo2', zlab='PCo3', 
            title='Trees embedded in dimension-reduced space',
            outfile=False,
        ):
        optioncheck(dimensions, [2,3])

        colours = 'bgrcmyk'
        colour_mapping = np.array([colours[i] for i in partition.partition_vector])
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

    from collection import Collection, Scorer
    from clustering import Clustering
    from lib.remote.utils import fileIO

    path_to_file = fileIO.path_to(__file__)
    test_data = fileIO.join_path(path_to_file, 'aa_alignments')

    c = Collection(
            input_dir=test_data,
            file_format='fasta',
            datatype='protein',
            compression='gz',
            )

    c.calc_NJ_trees()
    dm = c.distance_matrix('geo')

    cl = Clustering(dm)
    MDSdecomp = cl.MDS_decomp()
    spectralDecomp = cl.spectral_decomp(prune='estimate', local_scale='same')
    p = cl.MDS_cluster(4, spectralDecomp)
    print p

    coords2d = Plotter().decomp_to_coords(MDSdecomp, 2)
    coords3d = Plotter().decomp_to_coords(MDSdecomp, 3)

    fig2d = Plotter().embedding(coords2d, 2, p)
    fig3d = Plotter().embedding(coords3d, 3, p)
    plt.show()
