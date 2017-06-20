from __future__ import division
from builtins import object
import numpy as np
from scipy.spatial.distance import pdist, squareform

from .tasks import _fast_geo
from .constants import ISPY3

from tree_distance import PhyloTree

GOLDEN = (np.sqrt(5)-1)/2

# Functions for OutOfSampleMDS

def _eigen(mtx, inverse=False):
    if np.allclose(mtx, mtx.T):
        fn = np.linalg.eigh
    else:
        fn = np.linalg.eig
    vals, vecs = fn(mtx)
    ix = np.argsort(vals)[::-1]
    vals = vals[ix]
    vecs = vecs[:,ix]
    if inverse:
        return vecs, vals, np.linalg.inv(vecs)
    return vecs, vals


def _new_mean(mtx, sum_, size, new_row, index):
    old_row = mtx[index]
    return (sum_ - 2*old_row.sum() + old_row[index] + 2*new_row.sum() - new_row[index]) / size

def _new_rowmean(mtx, rowsum, nrow, new_row, index):
    old_row = mtx[index]
    tmp = rowsum - old_row + new_row
    tmp[index] = new_row.sum()
    return tmp / nrow

### Functions for OptimiseDistanceFit

def g_(x,a,c):
    """
    Calculate vector of equations of residuals,
    evaluated at x
    G[i] = Sum_j (x[j] - a[i,j])^2 - C[i]**2

    where:
    x is the point we're trying to fit (in M dimensions),
    a is the (N) already embedded points (in M dimensions),
    C is the expected distance
    """
    return ((x-a)**2).sum(1) - c**2

def g(x,a,c):
    """
    Christophe's suggestion for residuals,
    G[i] = Sqrt(Sum_j (x[j] - a[i,j])^2) - C[i] 
    """
    return np.sqrt(((x-a)**2).sum(1)) - c

def f(x, a, c):
    """ Objective function (sum of squared residuals) """
    v = g(x, a, c)
    return v.dot(v)

def jac_(x,a):
    """ Jacobian matrix of partial derivatives """
    return 2 * (x-a)

def jac(x,a):
    """ Jacobian matrix given Christophe's suggestion of f """
    return (x-a) / np.sqrt(((x-a)**2).sum(1))[:,np.newaxis]

def gradient(x, a, c):
    """ J'.G """
    return jac(x, a).T.dot(g(x, a, c))

def hessian(x, a):
    """ J'.J """
    j = jac(x, a)
    return j.T.dot(j)

def grad_desc_update(x, a, c, step=0.01):
    """ 
    Given a value of x, return a better x 
    using gradient descent 
    """
    return x - step * gradient(x,a,c)

def newton_update(x, a, c, step=1.0):
    """ 
    Given a value of x, return a better x 
    using newton-gauss 
    """
    return x - step*np.linalg.inv(hessian(x, a)).dot(gradient(x, a, c))

def levenberg_marquardt_update(x, a, c, damping=0.001):
    """ 
    Given a value of x, return a better x 
    using newton-gauss 
    """
    hess = hessian(x, a)
    return x - np.linalg.inv(hess + damping*np.diag(hess)).dot(gradient(x, a, c))

def golden_section_search(fn, a, b, tolerance=1e-5):
    """
    WIKIPEDIA IMPLEMENTATION
    golden section search
    to find the minimum of f on [a,b]
    f: a strictly unimodal function on [a,b]

    example:
    >>> f=lambda x:(x-2)**2
    >>> x=gss(f,1,5)
    >>> x
    2.000009644875678
    """
    c = b - GOLDEN*(b-a)
    d = a + GOLDEN*(b-a)
    while abs(c-d) > tolerance:
        fc, fd = fn(c), fn(d)
        if fc < fd:
            b = d
            d = c  #fd=fc;fc=f(c)
            c = b - GOLDEN*(b-a)
        else:
            a = c
            c = d  #fc=fd;fd=f(d)
            d = a + GOLDEN*(b-a)
    return (b+a)/2

def optimise_newton(x, a, c, tolerance=0.001):
    """
    Optimise value of x using newton gauss
    """
    x_new = x
    x_old = x-1 # dummy value
    while np.abs(x_new - x_old).sum() > tolerance:
        x_old = x_new
        x_new = newton_update(x_old, a, c)
    return x_new

def optimise_levenberg_marquardt(x, a, c, damping=0.001, tolerance=0.001):
    """
    Optimise value of x using levenberg-marquardt
    """
    x_new = x
    x_old = x-1 # dummy value
    f_old = f(x_new, a, c)
    while np.abs(x_new - x_old).sum() > tolerance:
        x_old = x_new
        x_tmp = levenberg_marquardt_update(x_old, a, c, damping)
        f_new = f(x_tmp, a, c)
        if f_new < f_old:
            damping = np.max(damping/10., 1e-20)
            x_new = x_tmp
            f_old = f_new
        else:
            damping *= 10.
    return x_new

def optimise_gradient_descent(x, a, c, tolerance=0.001):
    """
    Optimise value of x using gradient descent
    """
    x_new = x
    x_old = x-1 # dummy value
    while np.abs(x_new - x_old).sum() > tolerance:
        x_old = x_new
        step_size = golden_section_search(lambda step: f(grad_desc_update(x_old, a, c, step), a, c), -1.0, 1.0)
        x_new = grad_desc_update(x_old, a, c, step_size)
    return x_new

### Functions to add bootstraps to collections

def run_optimise_bootstrap_coords(boot_collection, ref_collection, ref_coords, task=_fast_geo, rooted=False, **kwargs):
    fit = np.empty((len(boot_collection), ref_coords.shape[1]))
    if ISPY3:
        query_trees = [PhyloTree(tree.encode(), rooted) for tree in boot_collection.trees]
        ref_trees = [PhyloTree(tree.encode(), rooted) for tree in ref_collection.trees]
    else:
        query_trees = [PhyloTree(tree, rooted) for tree in boot_collection.trees]
        ref_trees = [PhyloTree(tree, rooted) for tree in ref_collection.trees]
    for i, tree in enumerate(query_trees):
        ref_dists = np.array([task(tree, ref_tree, False) for ref_tree in ref_trees])
        opt = OptimiseDistanceFit(ref_coords.values, ref_dists)
        fit[i] = opt.newton(**kwargs)
    return fit

def run_out_of_sample_mds(boot_collection, ref_collection, ref_distance_matrix, index, dimensions, task=_fast_geo, rooted=False, **kwargs):
    """
    index = index of the locus the bootstrap sample corresponds to - only important if
            using recalc=True in kwargs
    """
    fit = np.empty((len(boot_collection), dimensions))
    if ISPY3:
        query_trees = [PhyloTree(tree.encode(), rooted) for tree in boot_collection.trees]
        ref_trees = [PhyloTree(tree.encode(), rooted) for tree in ref_collection.trees]
    else:
        query_trees = [PhyloTree(tree, rooted) for tree in boot_collection.trees]
        ref_trees = [PhyloTree(tree, rooted) for tree in ref_collection.trees]
    for i, tree in enumerate(query_trees):
        distvec = np.array([task(tree, ref_tree, False) for ref_tree in ref_trees])
        oos = OutOfSampleMDS(ref_distance_matrix)
        fit[i] = oos.fit(index, distvec, dimensions=dimensions, **kwargs)
    return fit

def run_analytical_fit(boot_collection, ref_collection, ref_coords, task=_fast_geo, rooted=False, **kwargs):
    fit = np.empty((len(boot_collection), ref_coords.shape[1]))
    if ISPY3:
        query_trees = [PhyloTree(tree.encode(), rooted) for tree in boot_collection.trees]
        ref_trees = [PhyloTree(tree.encode(), rooted) for tree in ref_collection.trees]
    else:
        query_trees = [PhyloTree(tree, rooted) for tree in boot_collection.trees]
        ref_trees = [PhyloTree(tree, rooted) for tree in ref_collection.trees]
    for i, tree in enumerate(query_trees):
        ref_dists = np.array([task(tree, ref_tree, False) for ref_tree in ref_trees])
        aft = AnalyticalFit(ref_coords.values, **kwargs)
        fit[i] = aft.fit(ref_dists)
    return fit

### Functions to assess closeness of fitted distances to reference

def stress(ref_cds, est_cds):
    """
    Kruskal's stress
    """
    ref_dists = pdist(ref_cds)
    est_dists = pdist(est_cds)
    return np.sqrt(((ref_dists - est_dists)**2).sum() / (ref_dists**2).sum())

def stress_dm(ref_distance_matrix, est_cds):
    ref_dists = squareform(ref_distance_matrix)
    est_dists = pdist(est_cds)
    return np.sqrt(((ref_dists - est_dists)**2).sum() / (ref_dists**2).sum())

def rmsd(ref_cds, est_cds):
    """
    Root-mean-squared-difference
    """
    ref_dists = pdist(ref_cds)
    est_dists = pdist(est_cds)
    return np.sqrt(((ref_dists - est_dists)**2).mean())

def rmsd_dm(ref_distance_matrix, est_cds):
    ref_dists = squareform(ref_distance_matrix)
    est_dists = pdist(est_cds)
    return np.sqrt(((ref_dists - est_dists)**2).mean())

### Classes


class OptimiseDistanceFit(object):
    """ 
    Usage: 
    opt = OptimiseDistanceFit(coords, dists, ['pairwise'|'adjacent']='adjacent'])
    opt.newton([start_x]=analytical_estimate, [tolerance]=1.0e-6)
    """
    def __init__(self, reference_coords, reference_dists, **kwargs):
        """
        Set up OptimiseDistanceFit with coordinates
        of reference points to fit to, and reference
        distances to find closest match to
        """
        self._a = reference_coords
        self._c = reference_dists
        self._analytical_fitter = AnalyticalFit(reference_coords, **kwargs)

    def residuals(self, x):
        """
        Calculate vector of equations of residuals evaluated
        at x
        """
        return g(x, self._a, self._c)

    def objective_fn(self, x):
        """
        Calculate objective function at x
        """
        return f(x, self._a, self._c)

    def jacobian(self, x):
        """
        Evaluate jacobian matrix at x
        """
        return jac(x, self._a)

    def gradient(self, x):
        """
        Evaluate gradient vector at x
        """
        return gradient(x, self._a, self._c)

    def hessian(self, x):
        """
        Evaluate hessian matrix at x
        """
        return hessian(x, self._a)

    def newton(self, start_x=None, tolerance=1.0e-6):
        """
        Optimise value of x using newton gauss
        """
        if start_x is None:
            start_x = self._analytical_fitter.fit(self._c)
        return optimise_newton(start_x, self._a, self._c, tolerance)

    def gradient_descent(self, start_x=None, tolerance=1.0e-6):
        """
        Optimise value of x using gradient descent
        """
        if start_x is None:
            start_x = self._analytical_fitter.fit(self._c)
        return optimise_gradient_descent(start_x, self._a, self._c, tolerance)

    def levenberg_marquardt(self, start_x=None, damping=1.0e-3, tolerance=1.0e-6):
        """
        Optimise value of x using levenberg marquardt
        """
        if start_x is None:
            start_x = self._analytical_fitter.fit(self._c)
        return optimise_levenberg_marquardt(start_x, self._a, self._c, tolerance)


class OutOfSampleMDS(object):
    def __init__(self, distance_matrix):
        """
        Store quantities calculated from a distance matrix,
        including the CMDS coordinate matrix
        """
        # Calculate all requirements once and store
        self.dmsq = distance_matrix**2
        self.rows, self.cols = distance_matrix.shape

        self.rowsum = (self.dmsq).sum(1)
        self.sum = (self.dmsq).sum()

        self.rowmean = (self.rowsum / self.rows)
        self.mean = self.sum / (self.rows*self.cols)

        B = -0.5 * (self.dmsq
                    - self.rowmean[np.newaxis]
                    - self.rowmean[:,np.newaxis]
                    + self.mean)
        
        U, l = _eigen(B)
        # l=l.astype(np.complex)
        l = np.clip(l, np.finfo(l.dtype).eps, np.inf)

        self.mult_factor = U*(1.0/np.sqrt(l))
        self.coords = U*np.sqrt(l)
        # self.mult_factor = self.mult_factor.astype(np.float)

    def new_B_row(self, row_index, distvec, recalc_mean=False):
        if recalc_mean:
            rowsum = distvec.sum()
            mean = _new_mean(self.dmsq, self.sum, self.rows*self.cols, distvec, row_index)
            rowmean = _new_rowmean(self.dmsq, self.rowsum, self.rows, distvec, row_index)
        else:
            mean = self.mean
            rowmean = self.rowmean
        b_row = -0.5 * (distvec - rowmean - rowmean[row_index] + mean)
        return b_row

    def new_coords(self, b_row):
        return np.dot(b_row, self.mult_factor)

    def fit(self, index, distvec, recalc=False, dimensions=3):
        """
        Replace distance matrix values at row/column index with
        distances in distvec, and compute new coordinates.
        Optionally use distvec to update means and (potentially) 
        get a better estimate.
        distvec values should be plain distances, not squared
        distances.
        """
        brow = self.new_B_row(index, distvec**2, recalc)
        return self.new_coords(brow)[:dimensions]


class AnalyticalFit(object):
    """
    Fit coords (x,y,[z]) so that distances from reference coordinates
    match closest to reference distances
    2D case - 
    reference coords = |a, b|
                       |c, d|
                       |e, f|
    reference dists  = [P, Q, R]
    Fit new coords (x, y):
     (x-a)^2 + (y-b)^2 = P^2   (eq.1) | This is an overdetermined
     (x-c)^2 + (y-d)^2 = Q^2   (eq.2) | system of simultaneous,
     (x-e)^2 + (y-f)^2 = R^2   (eq.3) | non-linear equations

    Construct overdetermined system of simultaneous *linear* equations
    by subtracting equations from each other - quadratic terms cancel
     (eq.4) = (eq.1) - (eq.2)
     (eq.5) = (eq.2) - (eq.3)
     (eq.6) = (eq.3) - (eq.1)
    
     (2c-2a)x + (2d-2b)y = P^2 - Q^2 - a^2 - b^2 + c^2 + d^2   (eq.4)
     (2e-2c)x + (2f-2d)y = Q^2 - R^2 - c^2 - d^2 + e^2 + f^2   (eq.5)
     (2a-2e)x + (2b-2f)y = R^2 - P^2 - e^2 - f^2 + a^2 + b^2   (eq.6)

    Best-fit solution found by least squares:
    (A is the matrix of coefficients on the left sides of eq.4--6,
     b is the vector of values on the rhs of eq.4--6
     A and the part of b that depends only on coords can be 
     calculated once and stored)

     Ax = b
     A'Ax = A'b
     x = (A'A)^-1 A'b
    """
    def __init__(self, ref_crds, method='adjacent'):
        """
        Construct A, part of b that depends on coords only,
        and Moore-Penrose pseudoinverse of A (i.e. (A'A)^-1A')
        """
        if method == 'adjacent':
            self._A, self._partial_b = self._make_A_and_part_of_b_adjacent(ref_crds)
            self._pinvA = np.linalg.pinv(self._A)
        elif method == 'pairwise':
            self._A, self._partial_b = self._make_A_and_part_of_b_pairwise(ref_crds)
            self._pinvA = np.linalg.pinv(self._A)
        else:
            raise ValueError('Unrecognised method {}'.format(method))
        self._method = method

    def _all_pairwise_comps(self, mtx):
        nrow = mtx.shape[0]
        ix_a, ix_b = np.triu_indices(nrow, 1)
        return mtx[ix_a], mtx[ix_b]

    def _rotate_rows(self, mtx):
        """
        rotate the matrix so all rows move up,
        and top row moves to bottom
        """
        return np.roll(mtx, -1, 0)

    def _make_A_and_part_of_b_adjacent(self, ref_crds):
        """
        Make A and part of b. See docstring of this class
        for answer to "What are A and b?"
        """
        rot = self._rotate_rows(ref_crds)
        A = 2*(rot - ref_crds)
        partial_b = (rot**2 - ref_crds**2).sum(1)
        return A, partial_b

    def _make_A_and_part_of_b_pairwise(self, ref_crds):
        m, n = self._all_pairwise_comps(ref_crds)
        A = 2*(n-m)
        partial_b = (n**2 - m**2).sum(1)
        return A, partial_b

    def _analytical_fit_adjacent(self, ref_dists):
        """
        Fit coords (x,y,[z]) so that distances from reference coordinates
        match closest to reference distances
        """
        dists = ref_dists**2
        rot_dists = self._rotate_rows(dists)
        b = dists - rot_dists + self._partial_b
        self._b = b
        return self._pinvA.dot(b)

    def _analytical_fit_pairwise(self, ref_dists):
        dists = ref_dists**2
        d1, d2 = self._all_pairwise_comps(dists)
        b = d1 - d2 + self._partial_b
        self._b = b
        return self._pinvA.dot(b)

    def fit(self, ref_dists):
        if self._method == 'adjacent':
            return self._analytical_fit_adjacent(ref_dists)
        elif self._method == 'pairwise':
            return self._analytical_fit_pairwise(ref_dists)
        else:
            raise ValueError('Unrecognised method {}'.format(method))


class BootstrapFitter(object):
    """
    Automate process of adding bootstraps to sample
    """
    pass
