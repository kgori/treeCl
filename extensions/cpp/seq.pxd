from  libcpp.string  cimport string as libcpp_string
from  libcpp cimport bool
cdef extern from "src/Seq.h" namespace "treeCl":
    cdef cppclass Seq:
        Seq() except +
        Seq(libcpp_string filename, libcpp_string file_format, libcpp_string datatype, bool interleaved) except +
        Seq(Seq& other) except +
        Seq operator+(Seq) except +
        void read_alignment(libcpp_string filename, libcpp_string file_format, libcpp_string datatype,
                            bool interleaved) except +
        void write_fasta(libcpp_string filename) except +
        void write_phylip(libcpp_string filename) except +
        size_t get_alignment_length() except +
        size_t get_number_of_sequences() except +
        Seq bootstrap_sample() except +
        void set_model(libcpp_string model_name) except +
        void set_gamma_rates(int ncat, double alpha) except +
        void compute_distances() except +
        void fast_compute_distances() except +
        void set_constant_rates(double value) except +
        bool is_dna() except +
        bool is_protein() except +
        libcpp_string bionj(bool fast) except +
        void optimize_numerical(double tolerance, long max_calls) except +
        void optimize_tree(double tolerance_before, double tolerance_during, long max_calls,
                           int nnis_per_round) except +
        void set_tree(libcpp_string newick) except +
        double get_likelihood() except +
        double get_likelihood(libcpp_string newick) except +
        void set_likelihood_object() except +
        libcpp_string newick() except +
        Seq simulate() except +
        # void operator+=(Seq) except +
