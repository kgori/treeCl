from __future__ import print_function

# standard library
import functools
import itertools

# third party
import scipy.spatial
from tree_distance import getEuclideanDistance, getGeodesicDistance, getRobinsonFouldsDistance,\
    getWeightedRobinsonFouldsDistance

# treeCl
from .utils import setup_progressbar

__all__ = ["eucdist", "eucdist_matrix", "geodist", "geodist_matrix", "rfdist", "rfdist_matrix", "wrfdist",
           "wrfdist_matrix"]


def _equalise_leaf_sets(t1, t2, inplace):
    intersect = t1 & t2
    if t1.labels != intersect:
        pruned1 = t1.prune_to_subset(intersect, inplace)
    else:
        pruned1 = t1
    if t2.labels != intersect:
        pruned2 = t2.prune_to_subset(intersect, inplace)
    else:
        pruned2 = t2
    return pruned1, pruned2


def _generic_distance_calc(fn, t1, t2, normalise, min_overlap=4):
    """(fn, t1, t2, normalise)

    Calculates the distance between trees t1 and t2. Can optionally be normalised to range [0, 1].
    If the trees have different leaf sets, the distance is calculated on their intersection. This incurs
    some overhead - if the trees are known to have the same leaves, then the underlying distance function [listed below]
    can be called instead, although in these cases the Tree arguments will need to be replaced with Tree.phylotree
    to make sure the appropriate data structure is passed.
    The distance is taken to be zero if the leaf overlap between the trees is less than `min_overlap`.
    E.g.
        treeCl.treedist.eucdist(t1, t2, False) is the leafset checking equivalent of
        treeCl.treedist.getEuclideanDistance(t1.phylotree, t2.phylotree, False)

    Distance functions:
        eucdist
        geodist
        rfdist
        wrfdist

    Underlying non-leafset-checking Distance functions:
        getEuclideanDistance
        getGeodesicDistance
        getRobinsonFouldsDistance
        getWeightedRobinsonFouldsDistance

    :param t1: Tree
    :param t2: Tree
    :param normalise: boolean
    :param min_overlap: int
    :return: float
    """
    if t1 ^ t2:
        if len(t1 & t2) < min_overlap:
            return 0
            #raise AttributeError('Can\'t calculate tree distances when tree overlap is less than two leaves')
        else:
            t1, t2 = _equalise_leaf_sets(t1, t2, False)

    return fn(t1.phylotree, t2.phylotree, normalise)


def _generic_matrix_calc(fn, trees, normalise, min_overlap=4):
    """(fn, trees, normalise)

    Calculates all pairwise distances between trees given in the parameter 'trees'.

    Distance functions:
        eucdist_matrix
        geodist_matrix
        rfdist_matrix
        wrfdist_matrix

    These wrap the leafset-checking functions. If the faster non-leafset-checking functions are needed, do this:
    scipy.spatial.distance(['getDistance'(t1.phylotree, t2.phylotree, normalise)
                                for (t1, t2) in itertools.combinations(trees, 2)])
    for your choice of 'getDistance' out of:
        getEuclideanDistance
        getGeodesicDistance
        getRobinsonFouldsDistance
        getWeightedRobinsonFouldsDistance

    :param trees: list or tuple, or some other iterable container type containing Tree objects
    :param normalise: boolean
    :param min_overlap: int
    :return: numpy.array
    """
    jobs = itertools.combinations(trees, 2)
    results = []
    pbar = setup_progressbar('Calculating tree distances', 0.5 * len(trees) * (len(trees) - 1))
    pbar.start()
    for i, (t1, t2) in enumerate(jobs):
        results.append(_generic_distance_calc(fn, t1, t2, normalise, min_overlap))
        pbar.update(i)
    pbar.finish()
    return scipy.spatial.distance.squareform(results)

eucdist = functools.partial(_generic_distance_calc, getEuclideanDistance)
geodist = functools.partial(_generic_distance_calc, getGeodesicDistance)
rfdist = functools.partial(_generic_distance_calc, getRobinsonFouldsDistance)
wrfdist = functools.partial(_generic_distance_calc, getWeightedRobinsonFouldsDistance)
eucdist_matrix = functools.partial(_generic_matrix_calc, getEuclideanDistance)
geodist_matrix = functools.partial(_generic_matrix_calc, getGeodesicDistance)
rfdist_matrix = functools.partial(_generic_matrix_calc, getRobinsonFouldsDistance)
wrfdist_matrix = functools.partial(_generic_matrix_calc, getWeightedRobinsonFouldsDistance)

eucdist.__doc__ = "eucdist" + _generic_distance_calc.__doc__
geodist.__doc__ = "geodist" + _generic_distance_calc.__doc__
rfdist.__doc__ = "rfdist" + _generic_distance_calc.__doc__
wrfdist.__doc__ = "wrfdist" + _generic_distance_calc.__doc__
eucdist_matrix.__doc__ = "eucdist_matrix" + _generic_matrix_calc.__doc__
geodist_matrix.__doc__ = "geodist_matrix" + _generic_matrix_calc.__doc__
rfdist_matrix.__doc__ = "rfdist_matrix" + _generic_matrix_calc.__doc__
wrfdist_matrix.__doc__ = "wrfdist_matrix" + _generic_matrix_calc.__doc__
