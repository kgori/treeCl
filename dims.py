#!/usr/bin/python

import numpy as np

"""
Some dimensional reduction functions for use on distance matricies
"""


def pcoordinate(nparray, *args):
    pass


def svd(nparray, *args):
    X = center_matrix(np.asarray(nparray))
    w1, U = np.linalg.eig(X.dot(X.T))
    U = U[:, w1.argsort()[::-1]]  # Sort vectors (columns) by eigenvalues (high to low)
    w2, V = np.linalg.eig(X.T.dot(X))
    V = V[:, w2.argsort()[::-1]]
    S = np.zeros((len(U), len(V)))
    np.fill_diagonal(S, sorted(np.sqrt(w1), reverse=True))
    return(U, S, V.T)


def center_matrix(nparray):
    rownorm = nparray - nparray.mean(axis=0)
    matrix = rownorm - nparray.mean(axis=1)[:, np.newaxis]
    return(matrix)
