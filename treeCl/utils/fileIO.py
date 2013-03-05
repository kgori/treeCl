#!/usr/bin/env python

__all__ = [
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
    'subprocess',
    'syscall',
    'verify',
    ]

import os
from subprocess import Popen, PIPE
import cPickle
import gzip

def basename(filename):
    return os.path.basename(filename)

def can_locate(filename):
    return (os.path.isfile(filename) if filename else False)


def can_open(directory):
    return (os.path.isdir(directory) if directory else False)


def delete(filename):
    if can_locate(filename):
        return os.remove(filename)
    else:
        raise FileError(filename)


def gpickle(obj, filename):

    if not filename.endswith('.gz'):
        filename += '.gz'
    cPickle.dump(obj, file=gzip.open(filename, 'wb'), protocol=-1)


def gunpickle(filename):

    return cPickle.load(gzip.open(filename, 'rb'))

def join_path(*elements):
    return os.path.join(*elements)


def locate_by_env(filename, path=None):
    path = os.getenv(path) or os.getenv('PATH', os.defpath)
    for directory in path.split(os.pathsep):
        if verify(filename, directory):
            print directory
            return os.path.abspath(directory)
        f = locate_by_dir(filename, directory)
        if f:
            return os.path.abspath(f)


def locate_by_dir(filename, directory=None):
    f = os.path.join(directory, filename)
    return (os.path.abspath(f) if can_locate(f) else None)


def locate_file(filename, env_var='', directory=''):
    f = locate_by_env(filename, env_var) or locate_by_dir(filename,
            directory)
    return (os.path.abspath(f) if can_locate(f) else None)


def path_to(filename):
    return os.path.dirname(filename)


def subprocess(cmd):
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    return process.communicate()


def syscall(cmd):
    return os.system(cmd)


def verify(filename, path):
    return can_locate(path) and os.path.basename(path) == filename
