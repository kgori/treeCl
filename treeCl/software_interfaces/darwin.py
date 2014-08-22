#!/usr/bin/env python
from __future__ import print_function

# standard library
from subprocess import Popen, PIPE

# third party
import numpy as np

# treeCl
from external import ExternalSoftware
from ..errors import filecheck, FileError, directorymake
from ..utils import fileIO


class Darwin(ExternalSoftware):
    """ Run commands through Darwin language"""

    default_binary = 'darwin'

    def __init__(self, tmpdir, verbosity=0):
        super(Darwin, self).__init__(tmpdir)
        self.tmpdir = directorymake(tmpdir)
        self.verbosity = verbosity
        self.comfile = '{0}/darcom.drw'.format(self.tmpdir)
        self.outfile = '{0}/output.drw'.format(self.tmpdir)

    def write(self, command):
        command = command.rstrip()
        writer = fileIO.fwriter(self.comfile)
        writer.write(command)
        if not command.endswith('\nquit;'):
            writer.write('\nquit;')
        writer.write('\n')
        writer.close()

    def execute(self, verbosity):
        try:
            from subprocess import DEVNULL  # py3k
        except ImportError:
            import os

            DEVNULL = open(os.devnull, 'wb')
        comfile = filecheck('{0}/darcom.drw'.format(self.tmpdir))
        p = Popen(self.binary, stdout=PIPE, stdin=PIPE, stderr=DEVNULL)
        cmd = "outf := \'{0}\'; ReadProgram(\'{1}\');".format(self.outfile,
                                                              comfile)
        if verbosity > 0:
            print(cmd)
        return p.communicate(input=cmd)[0]

    def read(self):
        filecheck(self.outfile)
        reader = fileIO.freader(self.outfile)
        result = reader.read()
        reader.close()
        return result

    def clean(self):
        for f in (self.comfile, self.outfile):
            try:
                fileIO.delete(f)
            except FileError:
                continue

    def run(self, command, verbosity=0):
        self.write(command)
        (stdout, stderr) = self.execute(verbosity)
        if verbosity > 0:
            print(stdout, stderr)
        result = self.read()
        self.clean()
        return result


def numpiser(s, dtype=None):
    elements = [line.strip().split() for line in s.strip().split('\n')]
    arr = np.array(elements, dtype=dtype)
    r, c = arr.shape
    if r == 1 or c == 1:
        arr = arr.reshape(max(r, c), )
    return arr


def runDarwin(cmd, tmpdir, dtype=None, **kwargs):
    dw = Darwin(tmpdir)
    output = dw.run(cmd, **kwargs)
    return numpiser(output, dtype)
