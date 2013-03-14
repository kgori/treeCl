#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

ext_modules = [Extension('evrot_extensions', ['evrot_extensions.c'])]

setup(name='evrot_extensions', version='0.0', ext_modules=ext_modules)
