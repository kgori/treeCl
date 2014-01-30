from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# from Cython.Build import cythonize
import pkg_resources

# data_dir = pkg_resources.resource_filename('autowrap', 'data_files')

extensions = Extension(name='tree_collection',
    sources=['py_wrapper.pyx',
             '../ProblemParser.cc',
             '../MinSqTree.cc',
             '../newick.cc',
             ],
    language='c++',
    # include_dirs=[data_dir],
    # extra_compile_args=['-std=c++11']
    )

setup(
    name = "tree_collection",
    cmdclass = { 'build_ext': build_ext},
    # ext_modules = [cythonize(extensions)],
    ext_modules = [extensions]
    )
