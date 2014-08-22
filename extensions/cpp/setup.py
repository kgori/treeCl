from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext

# from Cython.Build import cythonize
import pkg_resources

data_dir = pkg_resources.resource_filename('autowrap', 'data_files')

extensions = Extension(name='py_seq',
                       sources=['py_seq.pyx', 'src/Seq.cpp', 'src/SiteContainerBuilder.cpp', 'src/ModelFactory.cpp'],
                       language='c++',
                       include_dirs=[data_dir],
                       libraries=['bpp-core', 'bpp-seq', 'bpp-phyl'],
                       extra_compile_args=['-std=c++11'],
)

setup(
    name="py_seq",
    cmdclass={'build_ext': build_ext},
    # ext_modules = [cythonize(extensions)],
    ext_modules=[extensions]
)
