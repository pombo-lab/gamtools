"""

The resolution module
======================

This module contains functions for estimating the best resolution to use
for analysing a GAM dataset.

"""

import math

import numpy as np

# Some variable names are inherited from the SLICE library
# pylint: disable=invalid-name
from . import segregation


def segregation_info(segregation_table, skip_chroms):
    """
    Given a segregation table, retrieve the number of tubes,
    the diploid genome size and the bin size.
    """

    m = segregation_table.shape[1]

    segregation_windows = segregation_table.reset_index()[['chrom', 'start', 'stop']]
    segregation_windows['size'] = segregation_windows.stop - segregation_windows.start
    L = segregation_windows.loc[
        np.logical_not(segregation_windows.chrom.isin(skip_chroms)),
        'size'].sum() * 2

    bin_sizes = segregation_windows['size'].value_counts()
    frequent_sizes = [size for (size, freq)
                      in (bin_sizes / bin_sizes.sum()).iteritems()
                      if freq > 0.95]
    if not frequent_sizes:
        raise Exception('SLICE requires windows of even sizes')
    b = frequent_sizes[0]

    return m, L, b


def compute_heff_R(h, R, bL):

    return (h/R) + (2 * (bL ** (1/3)))


def get_detection_efficiency(segregation_table, h, R, skip_chroms=[], L=None, n_p=1, only_visible=True):
    """
    Given a GAM segregation table, calculate the efficiency of detection
    of individual loci. For example, 0.8 detection efficiency indicates that
    80% of the loci actually intersected by a slice are detected in our
    sequencing data.
    """

    m, genome_size, b = segregation_info(segregation_table, skip_chroms)

    if L is None:
        L = genome_size

    if only_visible:
        visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    else:
        visible = segregation_table

    m_visible = visible.shape[1]
    avg_m1s = visible.sum(axis=1).mean()
    avg_m1s_m = avg_m1s / m_visible

    coef = 1 / n_p
    bL = b / L

    heff_R = compute_heff_R(h, R, bL)

    return ((2 + heff_R) / heff_R) * (1 - (((1 - avg_m1s_m) ** coef) ** (1/2)))
