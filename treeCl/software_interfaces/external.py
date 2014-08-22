#!/usr/bin/env python
from __future__ import print_function
from ..errors import FileError
from ..utils import fileIO


class ExternalSoftware(object):
    """ Class for running external software. Wrappers for programs inherit from
    this class, and should implement their own read(), run() and write() methods
    """

    default_binary = ''
    default_env = ''
    local_dir = fileIO.path_to(__file__)

    def __init__(self, tmpdir, supplied_binary='', debug=False):

        self.flags = {}
        self.tempfiles = []

        # TRY TO FIND SUPPLIED BINARY
        path_to_supplied_binary = fileIO.locate_file(supplied_binary,
                                                     self.default_env, self.local_dir)

        if fileIO.can_locate(path_to_supplied_binary):
            self.binary = path_to_supplied_binary

        # FALL BACK TO DEFAULT BINARY
        else:
            default_binary = fileIO.locate_file(self.default_binary,
                                                self.default_env, self.local_dir)
            self.binary = default_binary

        if self.binary is None:
            raise FileError(supplied_binary)

        self.tmpdir = tmpdir.rstrip('/')
        self.debug = debug

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

    def call(self, verbose=False, dry_run=False):
        cmd = ' '.join([self.binary] + ['{0} {1}'.format(k, v) for (k, v) in
                                        self.flags.items()])
        if verbose:
            print(cmd)
        if dry_run:
            return cmd
        try:
            (stdout, stderr) = fileIO.subprocess(cmd)
        except KeyboardInterrupt:
            if not self.debug:
                self.clean()
                if verbose:
                    print('Cleaned up')
            raise
        if verbose:
            print(stdout, stderr)
        return (stdout, stderr)

    def clean(self):
        for fil in self.tempfiles:
            if fileIO.can_locate(fil):
                fileIO.delete(fil)


class TreeSoftware(ExternalSoftware):
    def __init__(self, record, tmpdir, supplied_binary=''):
        super(TreeSoftware, self).__init__(tmpdir, supplied_binary)
        self.record = record

    @property
    def record(self):
        return self._record

    @record.setter
    def record(self, sequence_record):
        self._record = sequence_record


class LSFJobHandlerException(Exception):
    pass


class LSFJobHandler(object):
    def __init__(self, job_handler, command_string, tmpdir):
        if (job_handler is None or command_string is None or tmpdir is None):
            raise LSFJobHandlerException('"None" argument given')

        self.job_handler = job_handler
        self.command_string = command_string
        self.tmpdir = tmpdir
        self.job_id = None

    def launch(self, job_name='lsf_job', verbose=False):
        self.job_id = bsub(job_name,
                           o='/dev/null',
                           e='/dev/null',
                           verbose=verbose)(self.command_string).job_id
        return self.job_id

    def clean(self):
        shutil.rmdir(self.tmpdir)
