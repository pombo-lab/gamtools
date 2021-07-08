from libcpp.string cimport string

cdef extern from "slice.h":
    int run_slice(string file_tube, string file_pi_out, string file_out_pi_thr,
                  string file_chr_names, string file_chr_indices,
                  int m, long L, int b, double h, double R, int n_p)

    long double compute_eps(long double avg_m1s_m, long double heff_R, int n_p)

    long double compute_heff_R(long double h, long double R, long double bL)

    int compilation_test()

def compute_detection_efficiency(avg_m1s_m, h, R, b, L, n_p):
    bL = b / L
    heff_R = compute_heff_R(h, R, bL)
    return compute_eps(avg_m1s_m, heff_R, n_p)

def slice(matrix_path,
          pi_out_path, threshold_out_path,
          chr_names_path, chr_indices_path,
          m, L, b, h, R, n_p):

    return run_slice(matrix_path.encode('UTF-8'),
                     pi_out_path.encode('UTF-8'),
                     threshold_out_path.encode('UTF-8'),
                     chr_names_path.encode('UTF-8'),
                     chr_indices_path.encode('UTF-8'),
                     m, L, b, h, R, n_p)

def test_compilation():
    return compilation_test()
