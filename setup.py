# -*- coding: utf-8 -*-
from __future__ import print_function

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, Extension

    def find_packages():
        return ['treeCl']
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

VERSION = '0.0.1'

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

extensions = [
    Extension(name='tree_collection',
              sources=[
                  'extensions/tree_collection/cython/py_wrapper.pyx',
                  'extensions/tree_collection/src/ProblemParser.cc',
                  'extensions/tree_collection/src/MinSqTree.cc',
                  'extensions/tree_collection/src/newick.cc',
              ],
              language='c++',
    ),
    # Extension(name='evrot_extensions',
    # sources=['extensions/evrot/evrot_extensions.c'],
    #     language='c',
    #     include_dirs=[numpy_get_include()]
    # ),
]

setup(name="treeCl",
      version=VERSION,
      author='Kevin Gori',
      author_email='kgori@ebi.ac.uk',
      description='Phylogenetic Clustering Package',
      url='https://github.com/kgori/treeCl.git',
      packages=find_packages(),
      include_package_data=True,
      package_data={
          'treeCl': ['software_interfaces/gtp.jar']
      },
      scripts=[
          'bin/simulator',
          'bin/treeCl',
          'bin/seqconvert',
          'bin/bootstrap',
      ],
      install_requires=[
          'biopython',
          'bsub',
          'cython',
          'dendropy',
          'fastcluster',
          'ggplot',
          'matplotlib',
          'numpy',
          'progressbar-latest',
          'scipy',
          'pandas',
          'scikit-learn',
      ],
      cmdclass={'build_ext': build_ext},
      ext_modules=extensions,
)
