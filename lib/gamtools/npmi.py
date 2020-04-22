"""
========================
The npmi module
========================

This module contains functions for calculating the normalized pointwise mutual
information (NPMI), which as of GAMtools v2.0 is the preferred metric for
generation of :ref:`proximity matrix <proximity_matrices>` as it does a better
job of normalizing for the differential detection of the two locations than our
previously recommended metric :func:`Dprime <get_dprime>`.


"""
# pylint: disable=invalid-name
import warnings

import numpy as np

def npmi_basic(A, B):
    """Calculate the normalized pointwise mutual information between A and B

    Takes two arrays (A and B) and returns their NPMI.

    This function is included mainly for debugging edge cases.
    """

    freq_A = sum(A)/len(A)
    freq_B = sum(B)/len(B)
    freq_AB = np.logical_and(A, B).sum()/len(A)
    if not freq_AB:
        return -1.0
    pmi = np.log(freq_AB/(freq_A * freq_B))
    npmi = pmi/-np.log(freq_AB)
    return npmi

def npmi_2d_slow(region1, region2):
    """Calculate an NPMI matrix between two regions

    Takes a list of :ref:`regions <regions>` and calculates the full NPMI
    matrix for all possible combinations of loci.

    This function is a pure python loop that calls
    :func:`npmi_basic <npmi_basic>`, so it is very slow but has the advantage
    of readability so is included for debugging edge cases.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the npmi \
            of all possible combinations of windows within the different regions.
    """

    npmi_matrix = np.zeros((len(region1), len(region2))) * np.NaN
    for i1, window_1_arr in enumerate(region1):
        for i2, window_2_arr in enumerate(region2):
            npmi_matrix[i1, i2] = npmi_basic(window_1_arr, window_2_arr)
    return npmi_matrix

def npmi_2d_fast(region1, region2):
    """Calculate an NPMI matrix between two regions

    Takes a list of :ref:`regions <regions>` and calculates the full NPMI
    matrix for all possible combinations of loci.

    This function uses vectorized numpy functions and is therefore much faster
    than npmi_2d_slow and is the recommended function for NPMI calculation.
    However, it is less readable and therefore potentially more difficult to
    debug.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the npmi \
            of all possible combinations of windows within the different regions.
    """

    M = region1.shape[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pXY = region1.dot(region2.T) / M
        pXpY = (region1.sum(1).reshape(-1, 1) / M) * (region2.sum(1) / M)
        pmi = np.log2(pXY / pXpY)
        npmi = pmi / -np.log2(pXY)
    return npmi

npmi_2d = npmi_2d_fast
