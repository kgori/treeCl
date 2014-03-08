try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, Extension
    def find_packages():
        return ['treeCl']
from Cython.Distutils import build_ext
from numpy import get_include as numpy_get_include
# from distutils.core import Extension
# from numpy.distutils.core import Extension as npExtension

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
    Extension(name='evrot_extensions',
        sources=['extensions/evrot/evrot_extensions.c'],
        language='c',
        include_dirs=[numpy_get_include()]
    ),
]


setup(name = "treeCl",
    version='1.0.0',
    author='Kevin Gori',
    author_email='kgori@ebi.ac.uk',
    description='TODO',
    url='https://github.com/kgori/gapmasker.git',
    packages=find_packages(),
    # packages=['treeCl'],
    scripts=[
        'bin/simulator',
        'bin/treeCl',
    ],
    install_requires=[
        'biopython',
        'bsub',
        'cython',
        'dendropy',
        'numpy',
        'scikit-learn',
    ],
    cmdclass = { 'build_ext': build_ext},
    ext_modules = extensions,
)

# setup(
#     name=name,
#     version='0.5.0',
#     author='Adrian Altenhoff',
#     author_email='adrian.altenhoff@inf.ethz.ch',
#     description='todoc',
#     packages=find_packages(),
#     install_requires=['lxml', 'progressbar-latest'],
#     scripts=['bin/familyanalyzer']
# )
