#!/usr/bin/env python

import os
import sys

__all__ = ['print_and_return', 'blank_line', 'restart_line', 'screen_width']


def print_and_return(s):
    blank_line()
    sys.stdout.write(str(s))
    sys.stdout.flush()
    restart_line()


def blank_line():
    sys.stdout.write(' ' * screen_width())
    sys.stdout.flush()
    restart_line()


def restart_line():
    sys.stdout.write('\r')
    sys.stdout.flush()


def screen_width():
    try:
        (_, c) = os.popen('stty size', 'r').read().split()
        return int(c)
    except:
        return 80
