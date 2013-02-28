#!/usr/bin/env python

from ..errors import FileError, filecheck_and_raise
from ..utils.fileIO import *

local_dir = path_to(__file__)


class ExternalSoftware(object):

    """ Class for running external software. Wrappers for programs inherit from
    this class, and should implement their own call(), read(), run() and write()
    methods """

    default_binary = ''
    default_env = ''

    def __init__(self, supplied_binary='', tmpdir='/tmp'):

        self.flags = {}
        self.tempfiles = []

        if can_locate(supplied_binary):
            self.binary = supplied_binary
        else:

            default_binary = locate_file(self.default_binary, self.default_env,
                    local_dir)
            self.binary = default_binary

        if self.binary is None:
            raise FileError(supplied_binary)

        self.tmpdir = tmpdir.rstrip('/')

    def __enter__(self):
        return self

    def __exit__(
        self,
        type,
        value,
        tb,
        ):
        self.clean()

    def __str__(self):
        desc = 'Wrapper for {0}'.format(self.default_binary)
        return desc

    def add_flag(self, flag, value):
        self.flags[flag] = value

    def add_tempfile(self, filename):
        self.tempfiles.append(filename)

    def remove_flag(self, flag):
        del self.flags[flag]

    def call(self, verbose=False):
        cmd = ' '.join([self.binary] + ['{0} {1}'.format(k, v) for (k, v) in
                       self.flags.items()])
        if verbose:
            print cmd
        (stdout, stderr) = subprocess(cmd)
        if verbose:
            print stdout, stderr
        return (stdout, stderr)

    def clean(self):
        for fil in self.tempfiles:
            if can_locate(fil):
                delete(fil)

    def read(self):
        pass

    def run(self):
        pass

    def write(self):
        pass
