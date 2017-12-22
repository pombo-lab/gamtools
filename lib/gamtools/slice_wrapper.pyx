from libcpp.string cimport string

cdef extern from "slice.h":
    int run_slice(string file_tube, string file_pi_out, string file_out_pi_thr, string file_chr_names, string file_chr_indices)
    int compilation_test()


def slice(segregation_file_path, pi_out_path, threshold_out_path, chr_names_path, chr_indices_path):

    return run_slice(segregation_file_path.encode('UTF-8'),
                     pi_out_path.encode('UTF-8'),
                     threshold_out_path.encode('UTF-8'),
                     chr_names_path.encode('UTF-8'),
                     chr_indices_path.encode('UTF-8'))

def test_compilation():
    return compilation_test()
