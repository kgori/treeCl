from  libcpp.string  cimport string as libcpp_string
from  libcpp.pair    cimport pair   as libcpp_pair
from  libcpp         cimport bool

cdef extern from "wrapper.cpp":
    cdef libcpp_pair[libcpp_string, double] compute(libcpp_string matrices,
                                                    libcpp_string mapping,
                                                    libcpp_string labels,
                                                    libcpp_string tree,
                                                    int iter,
                                                    bool loglik,
                                                    bool keep_topology,
                                                    bool quiet) except +
    cdef double fit(libcpp_string matrices,
                    libcpp_string mapping,
                    libcpp_string labels,
                    libcpp_string tree) except +
