from __future__ import print_function

# standard library
import itertools

# third party
import scipy

# treeCl
import tree

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


def _generic_distance_calc(fn, t1, t2, normalise):
    if t1 ^ t2:
        t1, t2 = _equalise_leaf_sets(t1, t2, False)
    if len(t1 & t2) < 2:
        raise AttributeError('Can\'t calculate tree distances when tree overlap is less than two leaves')
    return fn(t1.javatree, t2.javatree, normalise)


def _generic_matrix_calc(fn, trees, normalise):
    jobs = itertools.combinations(trees, 2)
    results = [_generic_distance_calc(fn, t1, t2, normalise) for (t1, t2) in jobs]
    return scipy.spatial.distance.squareform(results)


def eucdist(t1, t2, normalise=False):
    return _generic_distance_calc(tree.GTP.getEuclideanDistance, t1, t2, normalise)


def geodist(t1, t2, normalise=False):
    return _generic_distance_calc(tree.GTP.getGeodesicDistance, t1, t2, normalise)


def rfdist(t1, t2, normalise=False):
    return _generic_distance_calc(tree.GTP.getRobinsonFouldsDistance, t1, t2, normalise)


def wrfdist(t1, t2, normalise=False):
    return _generic_distance_calc(tree.GTP.getWeightedRobinsonFouldsDistance, t1, t2, normalise)


def eucdist_matrix(trees, normalise):
    return _generic_matrix_calc(tree.GTP.getEuclideanDistance, trees, normalise)


def geodist_matrix(trees, normalise):
    return _generic_matrix_calc(tree.GTP.getGeodesicDistance, trees, normalise)


def rfdist_matrix(trees, normalise):
    return _generic_matrix_calc(tree.GTP.getRobinsonFouldsDistance, trees, normalise)


def wrfdist_matrix(trees, normalise):
    return _generic_matrix_calc(tree.GTP.getWeightedRobinsonFouldsDistance, trees, normalise)
