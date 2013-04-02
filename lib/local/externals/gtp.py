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
            tmpf.write('\n'.join(tree.newick.rstrip() for tree in
                       trees))
        self.add_tempfile('{0}/geotrees.nwk'.format(self.tmpdir))
