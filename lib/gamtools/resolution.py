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
    genome_size = segregation_windows.loc[
        np.logical_not(segregation_windows.chrom.isin(skip_chroms)),
        'size'].sum() * 2

    bin_sizes = segregation_windows['size'].value_counts()
    frequent_sizes = [size for (size, freq)
                      in (bin_sizes / bin_sizes.sum()).iteritems()
                      if freq > 0.95]
    if not frequent_sizes:
        raise Exception('SLICE requires windows of even sizes')
    window_size = frequent_sizes[0]

    return m, genome_size, window_size


def compute_effective_slice_fraction(slice_thickness, nuclear_radius, window_fraction):

    return (slice_thickness/nuclear_radius) + (2 * (window_fraction ** (1/3)))


def get_detection_efficiency(segregation_table, slice_thickness, nuclear_radius, skip_chroms=[], genome_size=None, plexity=1, only_visible=True):
    """
    Given a GAM segregation table, calculate the efficiency of detection
    of individual loci. For example, 0.8 detection efficiency indicates that
    80% of the loci actually intersected by a slice are detected in our
    sequencing data.
    """

    m, table_genome_size, window_size = segregation_info(segregation_table, skip_chroms)

    if genome_size is None:
        genome_size = table_genome_size

    if only_visible:
        visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    else:
        visible = segregation_table

    m_visible = visible.shape[1]
    avg_m1s = visible.sum(axis=1).mean()
    avg_m1s_m = avg_m1s / m_visible

    coef = 1 / plexity
    window_fraction = window_size / genome_size

    effective_slice_fraction = compute_effective_slice_fraction(slice_thickness, nuclear_radius, window_fraction)

    return ((2 + effective_slice_fraction) / effective_slice_fraction) * (1 - (((1 - avg_m1s_m) ** coef) ** (1/2)))


def get_segregation_coverage_qc(segregation_table, skip_chroms=[],
                                only_visible=True, quantile=0.2):
    """
    Given a GAM segregation table, calculate x, such that 80% of the 
    windows in the segregation table are identified in an NP at least
    x times. For example, we recommend that in a "good quality" GAM
    dataset, 80% of the genomic windows will be detected at least 25
    times. 
    """

    if only_visible:
        visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    else:
        visible = segregation_table
        
    return visible.sum(axis=1).quantile(quantile)


def do_segregation_qc(segregation_table, slice_thickness, nuclear_radius, coverage_q=0.2,
                      skip_chroms=[], genome_size=None, plexity=1, only_visible=True):
    """
    Check that a GAM dataset has sufficient quality and depth.  Specifically,
    check that 80% of genomic windows have been detected at least 25 times
    and that the detection efficiency is at least 80%.
    """

    pass_qc = True

    window_coverage = get_segregation_coverage_qc(segregation_table, skip_chroms,
                                only_visible, coverage_q)

    if window_coverage >= 25:
        print('Coverage check passed. 80% of genomic windows are detected at '
              'least {:.0f} times'.format(window_coverage))
    else:
        print('Coverage check failed. 80% of genomic windows are detected at '
              'least {:.0f} times. This value should be at least 25, try increasing '
              'the genomic window size for this dataset.'.format(window_coverage))
        pass_qc = False

    detection_efficiency = get_detection_efficiency(segregation_table, slice_thickness, nuclear_radius,
                                                    skip_chroms, genome_size, plexity, only_visible)

    if detection_efficiency > 0.8:
        if detection_efficiency <= 1.02:
            print('Detection check passed. Genome-wide detection-efficiency '
                  '{:.1%}'.format(detection_efficiency))
        else:
            print('Detection check failed. Genome-wide detection-efficiency '
                  '{:.1%} - this cannot be higher than 100%, please check '
                  'your input parameters.'.format(
                   detection_efficiency))
            pass_qc = False
    else:
        print('Detection check failed. Genome-wide detection-efficiency '
              '{:.1%}. This value should be at least 80%.'.format(
               detection_efficiency))
        pass_qc = False

    return pass_qc

def resolution_check_from_args(args):
    """Extract parameters from an argparse namespace object and pass them
    to do_segregation_qc."""

    segregation_table = segregation.open_segregation(args.segregation_table)

    do_segregation_qc(segregation_table, args.slice_thickness, args.nuclear_radius,
                      skip_chroms=args.skip_chromosomes, genome_size=args.genome_size,
                      plexity=args.plexity, only_visible=args.only_visible)
