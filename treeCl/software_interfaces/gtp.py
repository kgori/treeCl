#!/usr/bin/env python

# standard library
import numbers
import os
import shutil
import tempfile
import time

# third party
from bsub import bsub
import numpy as np

# treeCl
from external import ExternalSoftware
from ..utils import fileIO


class LSFGTP(ExternalSoftware):

    default_binary = 'gtp.jar'
    default_env = 'GTP_PATH'
    local_dir = fileIO.path_to(__file__)

    def __init__(self, trees, tmpdir):
        super(LSFGTP, self).__init__(tmpdir)
        self.trees = trees
        self.temp_dirs = self.setup_temp_dirs()
        self.gtp_objects = self.setup_gtp_objects()

    @property
    def trees(self):
        return self._trees

    @trees.setter
    def trees(self, trees):
        self._trees = trees

    def write(self):
        pass

    def setup_temp_dirs(self):
        temp_dirs = [tempfile.mkdtemp(dir=self.tmpdir) for tree in self.trees]
        return temp_dirs

    def setup_gtp_objects(self):
        gtp_objects = [GTP(td) for td in self.temp_dirs]
        return gtp_objects

    def get_command_strings(self):
        return [self.gtp_objects[i].run(self.trees, i, dry_run=True)
                for i in range(len(self.trees))]

    def launch_lsf(self, command_strings, verbose=False, output='/dev/null'):
        curr_dir = os.getcwd()
        os.chdir(self.tmpdir)
        job_ids = [bsub('gtp_task',
                        o='/dev/null',
                        e='/dev/null',
                        verbose=verbose)(cmd).job_id
                   for cmd in command_strings]
        bsub.poll(job_ids)
        os.chdir(curr_dir)

    def read(self):
        rows = [gtp.read(len(self.trees), i)
                for i, gtp in enumerate(self.gtp_objects)]
        return np.vstack(rows)

    def clean(self):
        for gtp in self.gtp_objects:
            gtp.clean()
        for d in self.temp_dirs:
            shutil.rmtree(d)
        for f in os.listdir(self.tmpdir):
            if (f.endswith('.out') or f.endswith('.err')):
                os.remove(f)

    def run(self, verbose=False):
        command_strings = self.get_command_strings()
        self.launch_lsf(command_strings, verbose)
        time.sleep(2)
        matrix = self.read()
        return matrix

    def call(self):
        pass


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

    def check_row_value(self, row, maxval):
        if row is not None:
            if not isinstance(row, numbers.Number):
                raise ValueError('row parameter ({0}) is not a number'.format(
                    row))
            if row > (maxval - 1):
                raise ValueError('row value ({0}) is too big to index into'
                    ' a tree list of length {1}'.format(row, maxval))
            return int(row)

    def allrooted(self, trees):
        return all(tree.rooted for tree in trees)

    def call(self, rooted, row, dry_run=False):

        bincall = 'java -jar {0}'.format(self.binary)
        flags = []
        if row is not None:
            flags.append('-r')
        if not rooted:
            flags.append('-u')
        flags.append('-o')
        flags = ' '.join(flags)

        outf = '{0}/output.txt'.format(self.tmpdir)
        inf = '{0}/geotrees.nwk > /dev/null'.format(self.tmpdir)

        cmd = ' '.join((bincall, flags, outf, inf))
        if dry_run:
            return cmd
        fileIO.syscall(cmd)

    def pairwise(self, tree1, tree2):
        return self.run((tree1, tree2))[0, 1]

    def read(self, size, row, translation=None):
        row = self.check_row_value(row, size)
        translation = translation or self.translation
        if row is not None:
            matrix = np.zeros((1,size))
        else:
            matrix = np.zeros((size, size))

        try:
            with open('{0}/output.txt'.format(self.tmpdir)) as outf:
                for line in outf:
                    line = line.rstrip()
                    if line:
                        # print "line = ", line
                        (i, j, value) = line.split()
                        i = int(i)
                        j = translation[int(j)]
                        # print "i, j = ", i, j
                        value = float(value)
                        # print "value = ", value
                        # print "matrix dims = ", matrix.shape
                        matrix[i, j] = value
                        if row is None:
                            matrix[j, i] = value
            self.add_tempfile('{0}/output.txt'.format(self.tmpdir))
            return matrix
        except IOError, e:

            print 'There was an IOError: {0}'.format(e)
            print 'Geodesic distances couldn\'t be calculated'
            raise

    def run(self, trees, row=None, dry_run=False):
        row = self.check_row_value(row, len(trees))
        translation = self.write(trees, row)
        rooted = self.allrooted(trees)
        cmd = self.call(rooted, row, dry_run)
        if dry_run:
            return cmd
        try:
            matrix = self.read(len(trees), row, translation)
            self.clean()
            return matrix
        except IOError:
            print 'except'
            matrix = None
            raise

    def write(self, trees, row):
        row = self.check_row_value(row, len(trees))
        with open('{0}/geotrees.nwk'.format(self.tmpdir), 'w') as tmpf:
            if row is None:
                tmpf.write('\n'.join(self.get_safe_newick_string(tree)
                                     for tree in trees))
                translation = zip(range(len(trees)), range(len(trees)))

            else:
                tmpf.write(self.get_safe_newick_string(trees[row]) + '\n')
                translation = []
                for i, tree in enumerate(trees):
                    if i < row:
                        tmpf.write(self.get_safe_newick_string(tree) + '\n')
                        translation.append((i + 1, i))
                    elif i == row:
                        translation.append((0, i))
                    else:
                        tmpf.write(self.get_safe_newick_string(tree) + '\n')
                        translation.append((i, i))

        self.add_tempfile('{0}/geotrees.nwk'.format(self.tmpdir))
        self.translation = dict(translation)
        return self.translation # key:pos_new_list => value:pos_orig_list


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

def geodist(trees, tmpdir, row=None, lsf=False):

    if lsf and row is None:
        lsfgtp = LSFGTP(trees, tmpdir)
        return lsfgtp.run()
    else:
        g = GTP(tmpdir=tmpdir)
        l = len(trees)
        row = g.check_row_value(row, l)
        return g.run(trees, row)
