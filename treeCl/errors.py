#!/usr/bin/env python
from __future__ import print_function

# standard library
import numbers
import os
import sys

__all__ = [
    'FileError',
    'DirectoryError',
    'OptionError',
    'RangeError',
    'TreeBuildingError',
    'filecheck',
    'filequit',
    'directorymake',
    'directorycheck',
    'directoryquit',
    'optioncheck',
    'rangecheck',
    'isnumbercheck'
]


class FileError(Exception):
    """ Reporting when file errors occur """

    def __init__(self, filename=''):
        """ Store the file that gave the exception """

        self.value = filename

    def __str__(self):
        return 'Error opening file \'{0}\''.format(self.value)


class DirectoryError(FileError):
    """ Reporting when directory errors occur """

    def __str__(self):
        line1 = 'Error opening directory \'{0}\''.format(self.value)
        line2 = 'Directory doesn\'t exist'
        return '.'.join((line1, line2))


class OptionError(Exception):
    """ Reports on disallowed options passed to functions """

    def __init__(self, option, choices):
        self.value = option
        self.choices = choices

    def __str__(self):
        return '\'{0}\' is not a valid option. Valid options are {1}'.format(self.value,
                                                                             self.choices)


class RangeError(Exception):
    """ Raise exception when number is outside a valid range """

    def __init__(
            self,
            n,
            lower,
            upper,
    ):
        self.value = n
        self.lower = lower
        self.upper = upper

    def __str__(self):
        return '\'{0}\' is outside the valid range {1} - {2}'.format(self.value,
                                                                     self.lower, self.upper)

class TreeBuildingError(Exception):
    """ Raise exception when tree-building program fails """

    def __init__(self, error, program):
        self.error = error
        self.program = program

    def __str__(self):
        return '{0} failed: Error message:\n{1}'.format(self.program,
                                                        self.error)


def filecheck(filename):
    if not os.path.isfile(filename):
        raise FileError(filename)
    return filename


def filequit(filename):
    try:
        filecheck(filename)
    except FileError as e:
        print(e)
        sys.exit()


def directorycheck(directory):
    if not os.path.isdir(directory):
        raise DirectoryError(directory)
    return directory


def directoryquit(directory):
    try:
        if not os.path.isdir(directory):
            raise DirectoryError(directory)
    except DirectoryError as e:
        print(e)
        sys.exit()


def directorymake(directory, verbosity=0):
    try:
        directorycheck(directory)
    except DirectoryError as e:
        if verbosity > 1:
            print(e)
        if verbosity > 0:
            print('Creating \'{0}\''.format(directory))
        os.makedirs(directory)
    return directorycheck(directory)


def optioncheck(option, choices):
    if option not in choices:
        raise OptionError(option, choices)
    return option


def rangecheck(n, lower, upper):
    if n < lower or n > upper:
        raise RangeError(n, lower, upper)
    return n


def isnumbercheck(n):
    if not isinstance(n, numbers.Number):
        raise ValueError('{0} is not a number.')
    return n