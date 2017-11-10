#!/usr/bin/env python
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import range
from builtins import object
# standard library
import bz2
import pickle
import glob
import gzip
import io
import os
import shutil
import tempfile
from subprocess import Popen, PIPE

# treeCl
from ..errors import filecheck, directorycheck

__all__ = [
    'TempFile',
    'TempDir',
    'ChDir',
    'TempFileList'
    'basename',
    'can_locate',
    'can_open',
    'delete',
    'gpickle',
    'gunpickle',
    'join_path',
    'locate_by_dir',
    'locate_by_env',
    'locate_file',
    'path_to',
    'strip_extensions',
    'subprocess',
    'syscall',
    'verify',
]


class TempFile(object):

    def __init__(self, dir_=None, disable_delete=False):
        self.disable_delete = disable_delete
        if dir_ is not None and not os.path.exists(dir_):
            raise IOError('Directory "{}" does not exist'.format(dir_))
        self.dir = dir_

    def __enter__(self):
        self._fd, self._wrapped_tmp = tempfile.mkstemp(dir=self.dir)
        return os.path.abspath(self._wrapped_tmp)

    def __exit__(self, type, value, tb):
        # Always close the file, but only delete if disable_delete option allows it
        os.close(self._fd)
        if not self.disable_delete:
            os.remove(self._wrapped_tmp)


class TempDir(object):

    def __init__(self, disable_delete=False):
        self.disable_delete = disable_delete

    def __enter__(self):
        self._wrapped_tmpdir = tempfile.mkdtemp()
        return os.path.abspath(self._wrapped_tmpdir)

    def __exit__(self, type, value, tb):
        if not self.disable_delete:
            shutil.rmtree(self._wrapped_tmpdir)


class TempFileList(object):

    def __init__(self, filelist, disable_delete=False):
        self._filelist = filelist
        self.disable_delete = disable_delete

    def __enter__(self):
        return self._filelist

    def __exit__(self, type, value, tb):
        for fl in self._filelist:
            try:
                if not self.disable_delete:
                    os.remove(fl)
            except:
                pass  # No need to crash if deletion fails, just ignore


class ChDir(object):
    def __init__(self, working_dir):
        if not os.path.exists(working_dir):
            raise IOError('Directory "{}"" does not exist'.format(working_dir))
        self._cdir = os.getcwd()
        self._wdir = working_dir

    def __enter__(self):
        os.chdir(self._wdir)

    def __exit__(self, type, value, tb):
        os.chdir(self._cdir)


def basename(filename):
    return os.path.basename(filename)


def can_locate(filename):
    return os.path.isfile(filename) if filename else False


def can_open(directory):
    return os.path.isdir(directory) if directory else False


def delete(filename):
    return os.remove(filecheck(filename))


def delete_if_exists(filename):
    if os.path.isfile(filename):
        os.remove(filename)
        return True
    return False


def freader(filename, gz=False, bz=False):
    """ Returns a filereader object that can handle gzipped input """

    filecheck(filename)
    if filename.endswith('.gz'):
        gz = True
    elif filename.endswith('.bz2'):
        bz = True

    if gz:
        return gzip.open(filename, 'rb')
    elif bz:
        return bz2.BZ2File(filename, 'rb')
    else:
        return io.open(filename, 'rb')


def fwriter(filename, gz=False, bz=False):
    """ Returns a filewriter object that can write plain or gzipped output.
    If gzip or bzip2 compression is asked for then the usual filename extension will be added."""

    if filename.endswith('.gz'):
        gz = True
    elif filename.endswith('.bz2'):
        bz = True

    if gz:
        if not filename.endswith('.gz'):
            filename += '.gz'
        return gzip.open(filename, 'wb')
    elif bz:
        if not filename.endswith('.bz2'):
            filename += '.bz2'
        return bz2.BZ2File(filename, 'w')
    else:
        return open(filename, 'w')


def glob_by_extensions(directory, extensions):
    """ Returns files matched by all extensions in the extensions list """
    directorycheck(directory)
    files = []
    xt = files.extend
    for ex in extensions:
        xt(glob.glob('{0}/*.{1}'.format(directory, ex)))
    return files


def gpickle(obj, filename):
    if not filename.endswith('.gz'):
        filename += '.gz'
    pickle.dump(obj, file=gzip.open(filename, 'wb'), protocol=-1)


def gunpickle(filename):
    return pickle.load(gzip.open(filename, 'rb'))


def head(filename, n=10):
    """ prints the top `n` lines of a file """
    with freader(filename) as fr:
        for _ in range(n):
            print(fr.readline().strip())


def join_path(*elements):
    return os.path.join(*elements)


def locate_by_env(filename, path=None):
    path = os.getenv(path) or os.getenv('PATH', os.defpath)
    for directory in path.split(os.pathsep):
        if verify(filename, directory):
            print(directory)
            return os.path.abspath(directory)
        f = locate_by_dir(filename, directory)
        if f:
            return os.path.abspath(f)


def locate_by_dir(filename, directory=None):
    f = os.path.join(directory, filename)
    return os.path.abspath(f) if can_locate(f) else None


def locate_file(filename, env_var='', directory=''):
    """
    Locates a file given an environment variable or directory
    :param filename: filename to search for
    :param env_var: environment variable to look under
    :param directory: directory to look in
    :return: (string) absolute path to filename or None if not found
    """
    f = locate_by_env(filename, env_var) or locate_by_dir(filename, directory)
    return os.path.abspath(f) if can_locate(f) else None


def path_to(filename):
    return os.path.dirname(os.path.abspath(filename))


def strip_extensions(filename):
    toplevel = os.path.splitext(os.path.basename(filename))
    while toplevel[1] > '':
        toplevel = os.path.splitext(toplevel[0])
    return toplevel[0]


def subprocess(cmd):
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return process.communicate()


def syscall(cmd):
    return os.system(cmd)


def verify(filename, path):
    return can_locate(path) and os.path.basename(path) == filename
