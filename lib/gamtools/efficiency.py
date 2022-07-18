"""

The resolution module
======================

This module contains functions for estimating the best resolution to use
for analysing a GAM dataset.

"""

import math

# Some variable names are inherited from the SLICE library
# pylint: disable=invalid-name
from . import slice_wrapper, segregation, gam_slice


def get_detection_efficiency(segregation_file_path, h, R, skip_chroms=[], L=None, n_p=1):
    """
    Given a GAM segregation table, calculate the efficiency of detection
    of individual loci. For example, 0.8 detection efficiency indicates that
    80% of the loci actually intersected by a slice are detected in our
    sequencing data.
    """
    segregation_table = segregation.open_segregation(segregation_file_path)

    m, genome_size, b = gam_slice.segregation_info(segregation_table, skip_chroms)

    if L is None:
        L = genome_size

    visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    avg_m1s = visible.sum(axis=1).mean()
    avg_m1s_m = avg_m1s / m

    return slice_wrapper.compute_detection_efficiency(avg_m1s_m, h, R, b, L, n_p)


def compute_heff_R(h, R, bL):

    return (h/R) + (2 * (bL ** (1/3)))


def py_detection_efficiency(segregation_file_path, h, R, skip_chroms=[], L=None, n_p=1):
    """
    Given a GAM segregation table, calculate the efficiency of detection
    of individual loci. For example, 0.8 detection efficiency indicates that
    80% of the loci actually intersected by a slice are detected in our
    sequencing data.
    """
    segregation_table = segregation.open_segregation(segregation_file_path)

    m, genome_size, b = gam_slice.segregation_info(segregation_table, skip_chroms)

    if L is None:
        L = genome_size

    visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    avg_m1s = visible.sum(axis=1).mean()
    avg_m1s_m = avg_m1s / m

    coef = 1 / n_p
    bL = b / L
    heff_R = compute_heff_R(h, R, bL)

    return ((2 + heff_R) / heff_R) * (1 - (((1 - avg_m1s_m) ** coef) ** (1/2)))


#long double compute_eps(long double avg_m1s_m, const long double heff_R, int n_p ){
    #double coef=(double)(1.)/(double)(n_p);
    #return(  ((2+heff_R)/(heff_R)) * (1-sqrt(  pow((1-avg_m1s_m),coef)  ))  );
#}
