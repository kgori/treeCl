#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from builtins import str
from builtins import range
from builtins import object

# standard library
import itertools

# third party
import matplotlib.pyplot as plt
from matplotlib import cm as CM
from matplotlib.colors import hex2color
import numpy as np

# treeCl
from .distance_matrix import CoordinateMatrix, DistanceMatrix
from .partition import Partition
from .utils import flatten_list
from .colours import ggColorSlice

import logging
logger = logging.getLogger(__name__)


# Define some default sets of colours
SET2 = ["#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3"]
SET3 = ["#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f"]

def heatmap(dm, partition=None, cmap=CM.Blues, fontsize=10):
    """ heatmap(dm, partition=None, cmap=CM.Blues, fontsize=10)
    
    Produce a 2D plot of the distance matrix, with values encoded by
    coloured cells.

    Args:
        partition: treeCl.Partition object - if supplied, will reorder
                   rows and columns of the distance matrix to reflect
                   the groups defined by the partition
        cmap: matplotlib colourmap object  - the colour palette to use
        fontsize: int or None - sets the size of the locus lab

    Returns:
        matplotlib plottable object
    """
    assert isinstance(dm, DistanceMatrix)
    datamax = float(np.abs(dm.values).max())
    length = dm.shape[0]

    if partition:
        sorting = np.array(flatten_list(partition.get_membership()))
        new_dm = dm.reorder(dm.df.columns[sorting])
    else:
        new_dm = dm

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.xaxis.tick_top()
    ax.grid(False)

    tick_positions = np.array(list(range(length))) + 0.5
    if fontsize is not None:
        ax.set_yticks(tick_positions)
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(new_dm.df.columns, rotation=90, fontsize=fontsize, ha='center')
        ax.set_yticklabels(new_dm.df.index, fontsize=fontsize, va='center')

    cbar_ticks_at = [0, 0.5 * datamax, datamax]
    
    cax = ax.imshow(
        new_dm.values,
        interpolation='nearest',
        extent=[0., length, length, 0.],
        vmin=0,
        vmax=datamax,
        cmap=cmap,
    )
    cbar = fig.colorbar(cax, ticks=cbar_ticks_at, format='%1.2g')
    cbar.set_label('Distance')
    return fig


def _plotly_3d_scatter(coords, partition=None):
    """ _plotly_3d_scatter(coords, partition=None)

    Make a scatterplot of treeCl.CoordinateMatrix using the Plotly
    plotting engine
    """
    from plotly.graph_objs import Scatter3d, Data, Figure, Layout, Line, Margin, Marker
    # auto sign-in with credentials or use py.sign_in()

    colourmap = {
        'A':'#1f77b4', 
        'B':'#ff7f0e', 
        'C':'#2ca02c',
        'D':'#d62728',
        'E':'#9467bd',
        1:'#1f77b4', 
        2:'#ff7f0e', 
        3:'#2ca02c',
        4:'#d62728',
        5:'#9467bd'
    }

    df = coords.df
    if partition:
        assert len(partition.partition_vector) == df.shape[0]
        labels = [x+1 for x in partition.partition_vector]
    else:
        labels = [1 for _ in range(df.shape[0])]

    x, y, z = df.columns[:3]
    df['Label'] = labels
    
    colours = [colourmap[lab] for lab in df['Label']]
    trace = Scatter3d(x=df[x], y=df[y], z=df[z], mode='markers',
                      marker=Marker(size=9, color=colours, 
                                    line=Line(color=colours, width=0.5), opacity=0.8),
                      text=[str(ix) for ix in df.index])

    data = Data([trace])
    layout = Layout(
        margin=Margin(l=0, r=0, b=0, t=0 ),
        hovermode='x',
    )
    fig = Figure(data=data, layout=layout)
    return fig


def _add_sphere(ax):
    """ _add_sphere(ax)

    Add a wireframe unit sphere onto matplotlib 3D axes

    Args:
        ax - matplotlib 3D axes object

    Returns:
        updated matplotlib 3D axes
    """
    (u, v) = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = np.cos(u) * np.sin(v)
    y = np.sin(u) * np.sin(v)
    z = np.cos(v)
    ax.plot_wireframe(x, y, z, color='grey', linewidth=0.2)
    return ax


def plot_embedding(coordinates, partition=None, add_sphere=False, point_size=8,
        colours=None, labels=None, legend=True, outfile=False, **kwargs):
    """ plot_embedding(coordinates, partition=None, add_sphere=False, point_size=8,
    colours=None, labels=None, legend=True, outfile=False, **kwargs):
    
    Plot a 2D / 3D scatterplot of coordinates, optionally
    coloured by group membership.

    Args:
        coordinates: numpy array or treeCl.CoordinateMatrix -
            The coordinates of the points to plot. The number
            of columns determines the number of dimensions in
            the plot.
        add_sphere: bool -
            Add a wireframe sphere to a 3D plot. Spectral clustering
            places points on the surface of a unit sphere.
        colours: list of rgb hexes, or 'auto', or None -
            Colours to use to colour the points, as a list of
            RGB hex values. If None, defaults
            (colorbrewer set3). If 'auto', generates a set
            of colours equally spaced from the colour wheel.
        labels: Tuple(xlab, ylab, title, zlab) -
            Plot labels. Must be given in the above order.
            Missing options will be replaced by None. E.g.
            to set the title: (None, None, "Some points")
        outfile: str -
            Save figure to this filename
    """
    if isinstance(coordinates, CoordinateMatrix):
        coordinates = coordinates.values
    dimensions = min(3, coordinates.shape[1])
    partition = (partition or
                 Partition(tuple([0] * len(coordinates))))
    ngrp = partition.num_groups()

    if colours is None:
        colours = SET2
    elif colours == 'auto':
        colours = ggColorSlice(ngrp)

    colour_cycle = itertools.cycle(colours)
    colours = np.array([hex2color(c) for c in itertools.islice(colour_cycle, ngrp)])

    if labels is None:
        xlab, ylab, zlab, title = None, None, None, None
    else:
        if isinstance(labels, (tuple, list)):
            labels = list(labels[:4])
            labels.extend([None]*(4-len(labels)))
            xlab, ylab, title, zlab = labels

    fig = plt.figure()

    if dimensions == 3:
        ax = fig.add_subplot(111, projection='3d')
        if add_sphere:
            ax = _add_sphere(ax)
    else:
        ax = fig.add_subplot(111)
    
    members = partition.get_membership()
    for grp in range(ngrp):
        index = np.array(members[grp])
        points = coordinates[index,:dimensions].T
        ax.scatter(*points, s=point_size, c=colours[grp], edgecolor=colours[grp], label='Group {}'.format(grp+1), **kwargs)

    if xlab:
        ax.set_xlabel(xlab)
    if ylab:
        ax.set_ylabel(ylab)
    if zlab:
        ax.set_zlabel(zlab)
    if title:
        ax.set_title(title)
    if legend:
        plt.legend()
    if outfile:
        fig.savefig('{0}.pdf'.format(outfile))

    return fig


class Plotter(object):
    """ DEPRECATED
    """

    def __init__(self, *args, **kwargs):
        logger.warn("Plotter class is deprecated. Use module level functions\n"
                    "heatmap(...) and plot_embedding(...) instead.")

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

    def embedding_plotter(
            self, coordinates, partition=None, add_sphere=False, point_size=8,
            colours=None, labels=None, legend=True, outfile=False, **kwargs
    ):
        """ 
        Plot a 2D / 3D scatterplot of the coordinates, optionally
        coloured by group membership. 

        Parameters
        ==========
        coordinates [numpy array|treeCl.CoordinateMatrix] -
            The coordinates of the points to plot. The number
            of columns determines the number of dimensions in
            the plot.
        add_sphere [bool] -
            Add a wireframe sphere
        colours [None|list of rgb hexes|'auto'] -
            Colours to use to colour the points, as a list of
            RGB hex values. If None, defaults
            (colorbrewer set3). If 'auto', generates a set
            of colours similar to ggplot.
        labels [Tuple(xlab, ylab, title)] -
            Plot labels
        outfile [str] -
            Save figure
        """
        if isinstance(coordinates, CoordinateMatrix):
            coordinates = coordinates.values
        dimensions = min(3, coordinates.shape[1])
        partition = (partition or
                     Partition(tuple([0] * len(coordinates))))
        ngrp = partition.num_groups()

        if colours is None:
            colours = SET2
        elif colours == 'auto':
            colours = ggColorSlice(ngrp)

        colour_cycle = itertools.cycle(colours)
        colours = np.array([hex2color(c) for c in itertools.islice(colour_cycle, ngrp)])

        if labels is None:
            xlab, ylab, zlab, title = None, None, None, None
        else:
            if isinstance(labels, (tuple, list)):
                labels = list(labels[:4])
                labels.extend([None]*(4-len(labels)))
                xlab, ylab, zlab, title = labels

        fig = plt.figure()

        if dimensions == 3:
            ax = fig.add_subplot(111, projection='3d')
            if add_sphere:
                ax = self.sphere(ax)
        else:
            ax = fig.add_subplot(111)
        
        members = partition.get_membership()
        for grp in range(ngrp):
            index = np.array(members[grp])
            points = coordinates[index,:dimensions].T
            ax.scatter(*points, s=point_size, c=colours[grp], edgecolor=colours[grp], label='Group {}'.format(grp+1), **kwargs)

        if xlab:
            ax.set_xlabel(xlab)
        if ylab:
            ax.set_ylabel(ylab)
        if zlab:
            ax.set_zlabel(zlab)
        if title:
            ax.set_title(title)
        if legend:
            plt.legend()
        if outfile:
            fig.savefig('{0}.pdf'.format(outfile))

        return fig
