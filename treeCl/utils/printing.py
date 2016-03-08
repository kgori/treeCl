#!/usr/bin/env python
from __future__ import print_function
from builtins import str
from builtins import range
import sys

__all__ = ['print_and_return']


def print_and_return(s, stream=sys.stdout):
    print('\r\x1b[K{0}'.format(s), end='', file=stream)  # \x1b[K = Esc+[K = clear line
    stream.flush()


if __name__ == '__main__':
    # TEST
    import time

    for i in range(1, 11):
        print_and_return(i)
        time.sleep(0.3)
    print()

    for i in range(1, 21):
        print_and_return('.' * (22 - i - len(str(i))) + str(i))
        time.sleep(0.3)
    print()
