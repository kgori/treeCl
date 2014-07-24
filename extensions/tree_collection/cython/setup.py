from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext


extensions = Extension(name='tree_collection',
                       sources=['py_wrapper.pyx',
                                '../src/ProblemParser.cc',
                                '../src/MinSqTree.cc',
                                '../src/newick.cc',
                       ],
                       language='c++',
)

setup(
    name="tree_collection",
    cmdclass={'build_ext': build_ext},
    ext_modules=[extensions]
)
