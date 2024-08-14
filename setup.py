# -*- coding: utf-8 -*-
from __future__ import print_function

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, Extension

    def find_packages():
        return ['treeCl', 'treeCl.interfacing', 'treeCl.tasks', 'treeCl.utils']

from Cython.Distutils import build_ext
import pkg_resources
import platform
import re
import subprocess

# Facilities to install properly on Mac using clang
def is_clang(bin):
    proc = subprocess.Popen([bin, '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = str(b'\n'.join([stdout, stderr]).decode('ascii', 'ignore'))
    return not re.search(r'clang', output) is None

class my_build_ext(build_ext):
    def build_extensions(self):
        binary = self.compiler.compiler[0]
        if is_clang(binary):
            for e in self.extensions:
                e.extra_compile_args.append('-stdlib=libc++')

                if platform.system() == 'Darwin':
                    version_string, _, _ = platform.mac_ver()
                    version = tuple(int(n) for n in version_string.split('.'))

                    if version < (10, 9, 0):
                        # For very old Mac systems...
                        e.extra_compile_args.append('-mmacosx-version-min=10.7')
                        e.extra_link_args.append('-mmacosx-version-min=10.7')
                    else:
                        e.extra_compile_args.append('-mmacosx-version-min=10.9')
                        e.extra_compile_args.append('-stdlib=libc++')
                        e.extra_link_args.append('-mmacosx-version-min=10.9')
                        e.extra_link_args.append('-stdlib=libc++')

        build_ext.build_extensions(self)

compile_args = ['-std=c++11']

extensions = [
    Extension(name='tree_collection',
              sources=[
                  'extensions/tree_collection/cython/py_wrapper.pyx',
                  'extensions/tree_collection/src/ProblemParser.cc',
                  'extensions/tree_collection/src/MinSqTree.cc',
                  'extensions/tree_collection/src/newick.cc',
              ],
              language='c++',
              include_dirs=['extensions/tree_collection/src/eigen3'],
              extra_compile_args=compile_args,
    ),
]

# Install splash
VERSION = '0.1.41'

logo = """
═══════════ ╔═╗┬
┌┬┐┬─┐┌─┐┌─┐║  │
 │ ├┬┘├┤ ├┤ ╚═╝┴─┘
 ┴ ┴└─└─┘└─┘╭─────
┈┈┈┈┈┈┄┄┄┄┄─┤  ╭──
{versionfmt}╰──┤
══════════════ ╰──
""".format(versionfmt=VERSION.center(12))

print(logo)

setup(name="treeCl",
      version=VERSION,
      author='Kevin Gori',
      author_email='kcg25@cam.ac.uk',
      description='Phylogenetic Clustering Package',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      url='https://github.com/kgori/treeCl.git',
      packages=find_packages(),
      include_package_data=True,
      package_data={
          'treeCl': ['logging/logging.yaml']
      },
      scripts=[
          # 'bin/simulator',
          # 'bin/collapse',
          # 'bin/treeCl',
          # 'bin/seqconvert',
          # 'bin/bootstrap',
          # 'bin/npbs.py',
          # 'bin/pre_npbs.py',
      ],
      install_requires=[
          'biopython',
          'cython>=0.19.0',
          'dendropy>=4.0.0',
          'fastcluster',
          'futures; python_version == "2.7"',
          'ipython',
          'matplotlib',
          'nose',
          'numpy',
          'pandas',
          'phylo_utils==0.0.10',
          'progressbar-latest==2.4',
          'PyYaml',
          'scipy',
          'scikit-learn',
          'tree_distance>=1.0.12',
      ],
      cmdclass={'build_ext': my_build_ext},
      ext_modules=extensions,
      test_suite='tests',
)
