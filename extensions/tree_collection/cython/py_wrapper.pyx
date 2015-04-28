#cython: c_string_encoding=ascii  # for cython>=0.19
from  libcpp.string  cimport string as libcpp_string
# from  smart_ptr cimport shared_ptr
# from  AutowrapRefHolder cimport AutowrapRefHolder
from  libcpp cimport bool
from wrapper cimport compute as _compute_wrapper
from wrapper cimport fit as _fit_wrapper
# cdef extern from "autowrap_tools.hpp":
#     char * _cast_const_away(char *)

def compute(bytes matrices, bytes mapping, bytes labels, bytes tree, iter=5, loglik=True, keep_topology=False, quiet=True):
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
    :param keep_topology: If True, will assess the guidetree with no optimisation
    :type keep_topology: bool.
    :param quiet: If True, this will suppress printing progress to the screen
    :type quiet: bool.

    :returns: str -- the output tree, float -- the sum of squared residuals

    """
    assert isinstance(matrices, bytes), 'arg matrices wrong type'
    assert isinstance(mapping, bytes), 'arg mapping wrong type'
    assert isinstance(labels, bytes), 'arg labels wrong type'
    assert isinstance(tree, bytes), 'arg tree wrong type'
    assert isinstance(iter, (int, long)), 'arg iter wrong type'
    assert isinstance(loglik, (int, long)), 'arg loglik wrong type'
    assert isinstance(keep_topology, (int, long)), 'arg keep_topology wrong type'
    assert isinstance(quiet, (int, long)), 'arg quiet wrong type'

    _r = _compute_wrapper((<libcpp_string> matrices), (<libcpp_string> mapping), (<libcpp_string> labels),
                          (<libcpp_string> tree), (<int> iter), (<bool> loglik), (<bool> keep_topology), (<bool> quiet))
    cdef list py_result = [_r.first, _r.second]
    return py_result

def fit(bytes matrices, bytes mapping, bytes labels, bytes tree):
    assert isinstance(matrices, bytes), 'arg matrices wrong type'
    assert isinstance(mapping, bytes), 'arg mapping wrong type'
    assert isinstance(labels, bytes), 'arg labels wrong type'
    assert isinstance(tree, bytes), 'arg tree wrong type'

    cdef double _r = _fit_wrapper((<libcpp_string> matrices), (<libcpp_string> mapping), (<libcpp_string> labels),
                                  (<libcpp_string> tree))
    py_result = <double> _r
    return py_result
