# -*- coding: utf-8 -*-
from __future__ import print_function

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, Extension

    def find_packages():
        return ['treeCl', 'treeCl.interfacing', 'treeCl.tasks', 'treeCl.utils']
try:
    from Cython.Distutils import build_ext
except ImportError:
    print('You don\'t seem to have Cython installed.')
    print('Cython and Numpy are required for the')
    print('installation process, all other dependencies')
    print('will be installed automatically.')
    print('Install Cython [and Numpy] and try again.')
    import sys

    sys.exit()

try:
    from numpy import get_include as numpy_get_include
except ImportError:
    print('You don\'t seem to have Numpy installed.')
    print('Numpy and Cython are required for the')
    print('installation process, all other dependencies')
    print('will be installed automatically.')
    print('Install Numpy [and Cython] and try again.')
    import sys

    sys.exit()

import pkg_resources
import platform
import re
import subprocess

# Facilities to install properly on Mac using clang
def is_clang(bin):
    proc = subprocess.Popen([bin, '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = '\n'.join([stdout, stderr])
    return not re.search(r'clang', output) is None

class my_build_ext(build_ext):
    def build_extensions(self):
        binary = self.compiler.compiler[0]
        if is_clang(binary):
            for e in self.extensions:
                e.extra_compile_args.append('-stdlib=libc++')
                if platform.system() == 'Darwin':
                    e.extra_compile_args.append('-mmacosx-version-min=10.7')
        build_ext.build_extensions(self)

compile_args = ['-std=c++1y']

data_dir = pkg_resources.resource_filename("autowrap", "data_files")

extensions = [
    Extension(name='tree_collection',
              sources=[
                  'extensions/tree_collection/cython/py_wrapper.pyx',
                  'extensions/tree_collection/src/ProblemParser.cc',
                  'extensions/tree_collection/src/MinSqTree.cc',
                  'extensions/tree_collection/src/newick.cc',
              ],
              language='c++',
              include_dirs=[data_dir],
              extra_compile_args=['-std=c++11'],
    ),
]

# Install splash
VERSION = '0.1.0'

logo = """
═══════════ ╔═╗┬
┌┬┐┬─┐┌─┐┌─┐║  │
 │ ├┬┘├┤ ├┤ ╚═╝┴─┘
 ┴ ┴└─└─┘└─┘╭─────
┈┈┈┈┈┈┄┄┄┄┄─┤  ╭──
   V{0:s}   ╰──┤
══════════════ ╰──
""".format(VERSION)

print(logo)

setup(name="treeCl",
      version=VERSION,
      author='Kevin Gori',
      author_email='kgori@ebi.ac.uk',
      description='Phylogenetic Clustering Package',
      url='https://github.com/kgori/treeCl.git',
      packages=find_packages(),
      include_package_data=True,
      package_data={
          'treeCl': ['logging/logging.yaml']
      },
      scripts=[
          'bin/simulator',
          'bin/collapse',
          'bin/treeCl',
          'bin/seqconvert',
          'bin/bootstrap',
          'bin/npbs.py',
          'bin/pre_npbs.py',
      ],
      install_requires=[
          'autowrap',
          'biopython',
          'bpp',
          'cython',
          'dendropy',
          'fastcluster',
          'ipython',
          'matplotlib',
          'numpy',
          'pandas',
          'phylo_utils',
          'pllpy',
          'progressbar-latest',
          'scipy',
          'scikit-bio',
          'scikit-learn',
          'tree_distance',
      ],
      cmdclass={'build_ext': my_build_ext},
      ext_modules=extensions,
      test_suite='tests',
)
