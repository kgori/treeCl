#!/usr/bin/env python


import os
import sys
__all__ = [
    'FileError',
    'DirectoryError',
    'OptionError',
    'filecheck_and_raise',
    'filecheck_and_quit',
    'directorycheck_and_make',
    'directorycheck_and_raise',
    'directorycheck_and_quit',
    'optioncheck_and_raise',
    ]


class FileError(Exception):

    """
    Reporting when file errors occur
    """

    def __init__(self, filename=''):
        """
        Store the file that gave the exception
        """

        self.value = filename

    def __str__(self):
        return 'Error opening file \'{0}\''.format(self.value)


class DirectoryError(FileError):

    """
    Reporting when directory errors occur
    """

    def __str__(self):
        line1 = 'Error opening directory \'{0}\''.format(self.value)
        line2 = 'Directory doesn\'t exist'
        return '\n'.join((line1, line2))

class OptionError(Exception):

    """
    Reports on disallowed options passed to function
    """

    def __init__(self, option, choices):
        self.value = option
        self.choices = choices

    def __str__(self):
        return '{0} is not a valid option. Valid options are {1}'.format(self.value, self.choices)

def filecheck_and_raise(filename):
    if not os.path.isfile(filename):
        raise FileError(filename)
    return filename


def filecheck_and_quit(filename):
    try:
        filecheck_and_raise(filename)
    except FileError, e:
        print e
        sys.exit()


def directorycheck_and_raise(directory):
    if not os.path.isdir(directory):
        raise DirectoryError(directory)
    return directory


def directorycheck_and_quit(directory):
    try:
        if not os.path.isdir(directory):
            raise DirectoryError(directory)
    except DirectoryError, e:
        print e
        sys.exit()


def directorycheck_and_make(directory, verbose=True):
    try:
        directorycheck_and_raise(directory)
    except DirectoryError, e:
        if verbose:
            print e
            print 'Creating \'{0}\''.format(directory)
        os.makedirs(directory)
    return directorycheck_and_raise(directory)


def optioncheck_and_raise(option, choices):
    if option not in choices:
        raise OptionError(option, choices)
