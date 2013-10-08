#!/usr/bin/env python

from ...remote.externals.external import ExternalSoftware
from ...remote.errors import FileError, filecheck
from ...remote.utils import fileIO
import numpy as np


class GTP(ExternalSoftware):

    """
    Interacts with gtp.jar to calculate geodesic distances from trees
    """

    default_binary = 'gtp.jar'
    default_env = 'GTP_PATH'
    local_dir = fileIO.path_to(__file__)
    get_safe_newick_string = lambda self, tree: tree.as_string('newick', 
                internal_labels=False, suppress_rooting=True).rstrip()

    def __str__(self):
        desc = 'Wrapper for gtp.jar - geodesic distance calculator'
        authors = '(Owen, Megan, and J Scott Provan. 2011.'
        title = \
            'A Fast Algorithm for Computing Geodesic Distances in Tree Space,'
        doi = 'doi:10.1109/TCBB.2010.3)'

        details = \
            'Jarfile: {0}\nTemp directory: {1}'.format(self.binary,
                self.tmpdir)

        return '\n'.join((
            desc,
            authors,
            title,
            doi,
            '',
            details,
            '',
            ))

    def allrooted(self, trees):
        return all(tree.rooted for tree in trees)

    def call(self, rooted):

        bincall = 'java -jar {0}'.format(self.binary)
        if not rooted:
            flags = '-u -o'
        else:
            flags = '-o'
        outf = '{0}/output.txt'.format(self.tmpdir)
        inf = '{0}/geotrees.nwk > /dev/null'.format(self.tmpdir)

        cmd = ' '.join((bincall, flags, outf, inf))
        fileIO.syscall(cmd)

    def pairwise(self, tree1, tree2):
        return self.run((tree1, tree2))[0, 1]

    def read(self, size):
        matrix = np.zeros((size, size))

        try:
            with open('{0}/output.txt'.format(self.tmpdir)) as outf:
                for line in outf:
                    line = line.rstrip()
                    if line:
                        (i, j, value) = line.split()
                        i = int(i)
                        j = int(j)
                        value = float(value)
                        matrix[i, j] = matrix[j, i] = value
            self.add_tempfile('{0}/output.txt'.format(self.tmpdir))
            return matrix
        except IOError, e:

            print 'There was an IOError: {0}'.format(e)
            print 'Geodesic distances couldn\'t be calculated'
            raise

    def run(self, trees):
        self.write(trees)
        rooted = self.allrooted(trees)
        self.call(rooted)
        try:
            matrix = self.read(len(trees))
            self.clean()
            return matrix
        except IOError:
            print 'except'
            matrix = None
            raise

    def write(self, trees):
        with open('{0}/geotrees.nwk'.format(self.tmpdir), 'w') as tmpf:
            tmpf.write('\n'.join(self.get_safe_newick_string(tree) for tree in
                       trees))
        self.add_tempfile('{0}/geotrees.nwk'.format(self.tmpdir))


def breakpoints(n, m=100):
    """
    Takes an integer, returns k approx. equally sized integers,
    all < m, that sum to n
    """
    if n <= m:
        return n

    k = n / m + (1 if n % m > 0 else 0)
    s = n / k
    l = [s for _ in range(k)]

    r = n - sum(l)
    for i in range(r):
        l[i] += 1

    assert(sum(l) == n) # debug
    return np.array([0]+l).cumsum()

def offdiagonal_indices(breakpoints):
    matrix_length = breakpoints[-1]
    rows = []
    cols = []
    for i, curr in enumerate(breakpoints[1:], start=1):
        prev = breakpoints[i-1]
        for j in range(prev, curr):
            rows.extend([j]*(matrix_length-curr))
            cols.extend(range(curr,matrix_length))
    return rows, cols

def fill_diagonal_blocks(matrix, submatrices):
    i = 0
    j = 0
    for subm in submatrices:
        if subm.shape == (1,):
            r = c = 1
        else:
            r, c = subm.shape
        matrix[i:i+r, j:j+c] = subm
        i += r
        j += c
    return matrix

def geodist_diagonals(matrix, trees, breakpoints, calculator):
    """
    Calculator is GTP object
    """
    i = 0
    submatrices = []
    for b in breakpoints[1:]:
        tree_sublist = trees[i:b]
        i = b
        submatrices.append(calculator.run(tree_sublist))
    fill_diagonal_blocks(matrix, submatrices)

def geodist_offdiagonals(matrix, trees, breakpoints, calculator):
    for (r,c) in zip(*offdiagonal_indices(breakpoints)):
        matrix[r,c] = matrix[c,r] = calculator.pairwise(trees[r], trees[c])

def geodist(trees, tmpdir, breaks=None):
    g = GTP(tmpdir=tmpdir)
    if breaks is None:
        return g.run(trees)
    l = len(trees)
    bps = breakpoints(l, breaks)
    matrix = np.zeros((l,l))
    geodist_diagonals(matrix, trees, bps, g)
    geodist_offdiagonals(matrix, trees, bps, g)
    return matrix
