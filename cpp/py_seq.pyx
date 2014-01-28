#cython: c_string_encoding=ascii  # for cython>=0.19
from  libcpp.string  cimport string as libcpp_string
from  libcpp.set     cimport set as libcpp_set
from  libcpp.vector  cimport vector as libcpp_vector
from  libcpp.pair    cimport pair as libcpp_pair
from  libcpp.map     cimport map  as libcpp_map
from  smart_ptr cimport shared_ptr
from  AutowrapRefHolder cimport AutowrapRefHolder
from  libcpp cimport bool
from  libc.string cimport const_char
from cython.operator cimport dereference as deref, preincrement as inc, address as address
from seq cimport Seq as _Seq
cdef extern from "autowrap_tools.hpp":
    char * _cast_const_away(char *) 

cdef class Seq:

    cdef shared_ptr[_Seq] inst

    def __dealloc__(self):
         self.inst.reset()

    
    def set_constant_rates(self, double value ):
        assert isinstance(value, float), 'arg value wrong type'
    
        self.inst.get().set_constant_rates((<double>value))
    
    def bionj(self,  fast ):
        assert isinstance(fast, (int, long)), 'arg fast wrong type'
    
        cdef libcpp_string _r = self.inst.get().bionj((<bool>fast))
        py_result = <libcpp_string>_r
        return py_result
    
    def compute_distances(self):
        self.inst.get().compute_distances()
    
    def newick(self):
        cdef libcpp_string _r = self.inst.get().newick()
        py_result = <libcpp_string>_r
        return py_result
    
    def set_gamma_rates(self,  ncat , double alpha ):
        assert isinstance(ncat, (int, long)), 'arg ncat wrong type'
        assert isinstance(alpha, float), 'arg alpha wrong type'
    
    
        self.inst.get().set_gamma_rates((<int>ncat), (<double>alpha))
    
    def __add__(Seq self, Seq other not None):
        cdef _Seq  * this = self.inst.get()
        cdef _Seq * that = other.inst.get()
        cdef _Seq added = deref(this) + deref(that)
        cdef Seq result = Seq.__new__(Seq)
        result.inst = shared_ptr[_Seq](new _Seq(added))
        return result
    
    def write_fasta(self, str filename ):
        assert isinstance(filename, str), 'arg filename wrong type'
    
        self.inst.get().write_fasta((<libcpp_string>filename))
    
    def set_tree(self, str newick ):
        assert isinstance(newick, str), 'arg newick wrong type'
    
        self.inst.get().set_tree((<libcpp_string>newick))
    
    def bootstrap_sample(self):
        cdef _Seq * _r = new _Seq(self.inst.get().bootstrap_sample())
        cdef Seq py_result = Seq.__new__(Seq)
        py_result.inst = shared_ptr[_Seq](_r)
        return py_result
    
    def optimize_numerical(self, double tolerance ,  max_calls ):
        assert isinstance(tolerance, float), 'arg tolerance wrong type'
        assert isinstance(max_calls, (int, long)), 'arg max_calls wrong type'
    
    
        self.inst.get().optimize_numerical((<double>tolerance), (<long int>max_calls))
    
    def read_alignment(self, str filename , str file_format , str datatype ,  interleaved ):
        assert isinstance(filename, str), 'arg filename wrong type'
        assert isinstance(file_format, str), 'arg file_format wrong type'
        assert isinstance(datatype, str), 'arg datatype wrong type'
        assert isinstance(interleaved, (int, long)), 'arg interleaved wrong type'
    
    
    
    
        self.inst.get().read_alignment((<libcpp_string>filename), (<libcpp_string>file_format), (<libcpp_string>datatype), (<bool>interleaved))
    
    def optimize_tree(self, double tolerance_before , double tolerance_during ,  max_calls ,  nnis_per_round ):
        assert isinstance(tolerance_before, float), 'arg tolerance_before wrong type'
        assert isinstance(tolerance_during, float), 'arg tolerance_during wrong type'
        assert isinstance(max_calls, (int, long)), 'arg max_calls wrong type'
        assert isinstance(nnis_per_round, (int, long)), 'arg nnis_per_round wrong type'
    
    
    
    
        self.inst.get().optimize_tree((<double>tolerance_before), (<double>tolerance_during), (<long int>max_calls), (<int>nnis_per_round))
    
    def get_alignment_length(self):
        cdef size_t _r = self.inst.get().get_alignment_length()
        py_result = <size_t>_r
        return py_result
    
    def is_protein(self):
        cdef bool _r = self.inst.get().is_protein()
        py_result = <bool>_r
        return py_result
    
    def is_dna(self):
        cdef bool _r = self.inst.get().is_dna()
        py_result = <bool>_r
        return py_result
    
    def simulate(self):
        cdef _Seq * _r = new _Seq(self.inst.get().simulate())
        cdef Seq py_result = Seq.__new__(Seq)
        py_result.inst = shared_ptr[_Seq](_r)
        return py_result
    
    def write_phylip(self, str filename ):
        assert isinstance(filename, str), 'arg filename wrong type'
    
        self.inst.get().write_phylip((<libcpp_string>filename))
    
    def fast_compute_distances(self):
        self.inst.get().fast_compute_distances()
    
    def __copy__(self):
       cdef Seq rv = Seq.__new__(Seq)
       rv.inst = shared_ptr[_Seq](new _Seq(deref(self.inst.get())))
       return rv
    
    def _init_0(self):
        self.inst = shared_ptr[_Seq](new _Seq())
    
    def _init_1(self, str filename , str file_format , str datatype ,  interleaved ):
        assert isinstance(filename, str), 'arg filename wrong type'
        assert isinstance(file_format, str), 'arg file_format wrong type'
        assert isinstance(datatype, str), 'arg datatype wrong type'
        assert isinstance(interleaved, (int, long)), 'arg interleaved wrong type'
    
    
    
    
        self.inst = shared_ptr[_Seq](new _Seq((<libcpp_string>filename), (<libcpp_string>file_format), (<libcpp_string>datatype), (<bool>interleaved)))
    
    def __init__(self, *args):
        if not args:
             self._init_0(*args)
        elif (len(args)==4) and (isinstance(args[0], str)) and (isinstance(args[1], str)) and (isinstance(args[2], str)) and (isinstance(args[3], (int, long))):
             self._init_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def get_number_of_sequences(self):
        cdef size_t _r = self.inst.get().get_number_of_sequences()
        py_result = <size_t>_r
        return py_result
    
    def set_model(self, str model_name ):
        assert isinstance(model_name, str), 'arg model_name wrong type'
    
        self.inst.get().set_model((<libcpp_string>model_name))
    
    def _get_likelihood_0(self):
        cdef double _r = self.inst.get().get_likelihood()
        py_result = <double>_r
        return py_result
    
    def _get_likelihood_1(self, str newick ):
        assert isinstance(newick, str), 'arg newick wrong type'
    
        cdef double _r = self.inst.get().get_likelihood((<libcpp_string>newick))
        py_result = <double>_r
        return py_result
    
    def get_likelihood(self, *args):
        if not args:
            return self._get_likelihood_0(*args)
        elif (len(args)==1) and (isinstance(args[0], str)):
            return self._get_likelihood_1(*args)
        else:
               raise Exception('can not handle type of %s' % (args,))
    
    def set_likelihood_object(self):
        self.inst.get().set_likelihood_object() 

