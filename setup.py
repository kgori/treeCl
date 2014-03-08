# from distutils.core import setup
# from distutils.extension import Extension
from Cython.Distutils import build_ext
from numpy.distutils.core import setup, Extension

extensions = [
    Extension(name='tree_collection',
        sources=['py_wrapper.pyx',
            '../src/ProblemParser.cc',
            '../src/MinSqTree.cc',
            '../src/newick.cc',
        ],
        language='c++',
    ),
    Extension(name='evrot_extensions',
        sources=['evrot_extensions.c'],
        language='c',
    ),
]


setup(name = "tree_collection",
    version='1.0.0',
    author='Kevin Gori',
    author_email='kgori@ebi.ac.uk',
    description='TODO',
    url='https://github.com/kgori/gapmasker.git',
    packages=['treeCl',
    ],
    scripts=['bin/simulator',
        'bin/treeCl',
    ],
    install_requires=['dendropy', 'numpy', 'scikit-learn', 'biopython', 'bsub'],
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
