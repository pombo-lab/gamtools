#ifndef _functions_h
#define _functions_h

#include <iostream>
using namespace std;

int run_slice(string file_tube, string file_pi_out, string file_out_pi_thr,
              string file_chr_names, string file_chr_indices,
              int m, long L, int b, double h, double R, int n_p);

long double compute_eps(long double avg_m1s_m, long double heff_R, int n_p);

long double compute_heff_R(long double h, long double R, long double bL);

int compilation_test();

#endif
