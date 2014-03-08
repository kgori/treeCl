#cython: c_string_encoding=ascii  # for cython>=0.19
from  libcpp.string  cimport string as libcpp_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from wrapper cimport compute as _compute_wrapper

def compute(str matrices , str mapping , str labels , str tree ,  iter=3 ,  quiet=True ):
    """ This is a wrapper to the MinSquareTreeCollection::compute function

    :param matrices: Matrix of distances (upper triangle) and variances (lower triangle)
    :type matrices: str.
    :param mapping: Genome map
    :type mapping: str.
    :param labels: A list of species labels
    :type labels: str.
    :param tree: The guide tree in newick format
    :type tree: str.
    :param iter: Number of iterations to do in the compute (resets on each tree swap)
    :type iter: int.
    :param quiet: If True, this will suppress printing progress to the screen
    :type quiet: bool.

    :returns: str -- the output tree, float -- the sum of squared residuals

    """
    assert isinstance(matrices, str), 'arg matrices wrong type'
    assert isinstance(mapping, str), 'arg mapping wrong type'
    assert isinstance(labels, str), 'arg labels wrong type'
    assert isinstance(tree, str), 'arg tree wrong type'
    assert isinstance(iter, (int, long)), 'arg iter wrong type'
    assert isinstance(quiet, (int, long)), 'arg quiet wrong type'

    _r = _compute_wrapper((<libcpp_string>matrices), (<libcpp_string>mapping), (<libcpp_string>labels), (<libcpp_string>tree), (<int>iter), (<bool>quiet))
    cdef tuple py_result = (_r.first, _r.second)
    return py_result
