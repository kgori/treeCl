#!/usr/bin/env python

import numpy as np

"""
Some dimensional reduction functions for use on distance matricies
"""


def pcoordinate(nparray, *args):
    pass


def svd(nparray):
    """
    Use np.linalg.svd - this doesn't have same sign-ambiguity problem
    as this method (+ and - eigenvectors both satisfy characteristic eq.)
    """
    X = double_centre(np.asarray(nparray))
    X = nparray
    w1, U = np.linalg.eig(X.dot(X.T))
    U = U[:, w1.argsort()[::-1]]  # Sort vectors (columns) by eigenvalues (high to low)
    w2, V = np.linalg.eig(X.T.dot(X))
    V = V[:, w2.argsort()[::-1]]
    S = np.zeros((len(U), len(V)))
    np.fill_diagonal(S, sorted(np.sqrt(w1), reverse=True))
    return(U, S, V.T)

def gram_schmidt(a):
    """
    Calculates an orthonormalised matrix of (column) vectors from
    input matrix a.
    Orthonormal vectors have length = 1 and dot product = 0 (perpendicular)
    """
    nrows, ncols = a.shape

    # set up empty matrix to hold orthonormalised vecs
    u = np.zeros((nrows, ncols))

    # fill first column
    u[:, 0] = a[:, 0] / np.linalg.norm(a[:, 0])

    # k-th orthogonal vector given by
    # $w_k = v_k - \sum_{i=1}^{k-1} u_i \cdot v_k \times u_i$
    # k-th orthonormal vector is w_k / ||w_k||
    for k in range(1, ncols):
        v_k = a[:, k]
        sum_ = 0
        for i in range(k):
            u_i = u[:, i]
            sum_ += u_i.dot(v_k) * u_i
        w_k = v_k - sum_
        w_k /= np.linalg.norm(w_k)
        u[:, k] = w_k

    return u

def double_centre(D, square_input=True):
    """
    Double-centres the input matrix:
    From each element:
    Subtract the row mean
    Subtract the column mean
    Add the grand mean
    Divide by -2
    Method from: Torgerson, W S (1952). Multidimensional scaling: I. Theory and
    method.
    Alternatively M = -0.5 * (I - 1/n)D[^2](I - 1/n) (slow for large D)
    """

    M = (D*D if square_input else D)
    (rows, cols) = M.shape

    cm = np.mean(M, axis=0)  # column means
    rm = np.mean(M, axis=1).reshape((rows, 1))  # row means
    gm = np.mean(cm)  # grand mean
    M -= rm + cm - gm
    M /= -2
    return M


def center_matrix(nparray):
    rownorm = nparray - nparray.mean(axis=0)
    matrix = rownorm - nparray.mean(axis=1)[:, np.newaxis]
    return(matrix)
