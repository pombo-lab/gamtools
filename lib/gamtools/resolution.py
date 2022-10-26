"""

The resolution module
======================

This module contains functions for estimating the best resolution to use
for analysing a GAM dataset.

Two metrics are used:
 * Detection efficiency
 * Window coverage

Detection efficiency
--------------------

Detection efficiency is an estimate of the percentage of DNA actually present
in each library that was detected after WGA and sequencing. For example if your
sample contains three NPs with a total of 10 million base pairs of DNA, but the
detected positive windows only cover 6 million base pairs of DNA, the detection
efficiency would be 60%. This value is estimated for a whole dataset and we
would recommend it should be 70% or higher. If detection efficiency is too low,
this can be a sign that your window size is too small (i.e. you are trying to
work at too high a resolution). Datasets with lower detection efficiency can
still be usable, but a greater number of libraries will need to be sequenced.

Window coverage
---------------

In order to have good statistical power, each genomic window should be observed
in the dataset with a sufficient frequency. We recommend that 80% of the windows
in the genome should be detected at least 20 times. If this number is too low,
increasing the window size may help a little but you may need to collect more
libraries.

"""

import numpy as np

from . import segregation


def segregation_info(segregation_table, skip_chroms):
    """
    Given a segregation table, retrieve the number of tubes,
    the diploid genome size and the bin size.

    :param segregation_table: Input :ref:`segregation table <segregation_table>`
    :param list skip_chroms: List of chromosome names to exclude.
    :returns: Number of libraries, total size of the genome and size of each window
    """

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

    return genome_size, window_size


def compute_eff_slice_fraction(slice_thickness, nuclear_radius, window_fraction):
    """
    Compute the ratio between slice thickness and nuclear radius, taking into account
    the window size. This is important because a genomic window does not have to fall
    entirely within the NP to be detected. The larger the window size, the further DNA
    can fall outside of the physical section and still be 'effectively' within the NP.

    :param float slice_thickness: Thickness of the NP.
    :param float nuclear_radius: Average nuclear radius.
    :param float window_fraction: Genomic window size divided by total genome size.
    :returns: Effective slice thickness
    """

    return (slice_thickness/nuclear_radius) + (2 * (window_fraction ** (1/3)))


def get_detection_efficiency(segregation_table, slice_thickness, #pylint: disable=too-many-arguments
                             nuclear_radius, skip_chroms=None, genome_size=None,
                             plexity=1, only_visible=True):
    """
    Given a GAM segregation table, calculate the efficiency of detection
    of individual loci. For example, 0.8 detection efficiency indicates that
    80% of the loci actually intersected by a slice are detected in our
    sequencing data.

    :param segregation_table: Input :ref:`segregation table <segregation_table>`
    :param float slice_thickness: Thickness of the NP.
    :param float nuclear_radius: Average nuclear radius.
    :param list skip_chroms: List of chromosome names to exclude.
    :param float genome_size: Total size of the genome (if None, calculate from \
          the segregation table.
    :param int plexity: Number of nuclear profiles in each sample.
    :param bool only_visible: If true, don't include genomic windows that are \
          never detected across the dataset in calculation. This is the default \
          behaviour.
    :returns: Number of libraries, total size of the genome and size of each window
    """

    if skip_chroms is None:
        skip_chroms = []

    table_genome_size, window_size = segregation_info(segregation_table, skip_chroms)

    if genome_size is None:
        genome_size = table_genome_size

    if only_visible:
        visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    else:
        visible = segregation_table

    avg_m1s_m =  visible.sum(axis=1).mean() / visible.shape[1]

    coef = 1 / plexity
    window_fraction = window_size / genome_size

    effective_slice_fraction = compute_eff_slice_fraction(
        slice_thickness, nuclear_radius, window_fraction)

    return ((2 + effective_slice_fraction) / effective_slice_fraction) * (1 -
            (((1 - avg_m1s_m) ** coef) ** (1/2)))


def get_segregation_coverage_qc(segregation_table, skip_chroms=None,
                                only_visible=True, quantile=0.2):
    """
    Given a GAM segregation table, calculate x, such that 80% of the
    windows in the segregation table are identified in an NP at least
    x times. For example, we recommend that in a "good quality" GAM
    dataset, 80% of the genomic windows will be detected at least 20
    times.

    :param segregation_table: Input :ref:`segregation table <segregation_table>`
    :param list skip_chroms: List of chromosome names to exclude.
    :param bool only_visible: If true, don't include genomic windows that are \
          never detected across the dataset in calculation. This is the default \
          behaviour.
    :param float quantile: The quantile at which to calculate window coverage. \
          The default value of 0.2 returns the coverage exceeded by 80% of \
          genomic windows.
    :returns: Number of times 80% of genomic windows are detected.
    """

    if skip_chroms is None:
        skip_chroms = []

    if only_visible:
        visible = segregation_table.loc[segregation_table.sum(axis=1) > 0]
    else:
        visible = segregation_table

    return visible.sum(axis=1).quantile(quantile)


def do_segregation_qc(segregation_table, slice_thickness, nuclear_radius, coverage_q=0.2, #pylint: disable=too-many-arguments
                      skip_chroms=None, genome_size=None, plexity=1, only_visible=True):
    """
    Check that a GAM dataset has sufficient quality and depth.  Specifically,
    check that 80% of genomic windows have been detected at least 20 times
    and that the detection efficiency is at least 70%.

    :param segregation_table: Input :ref:`segregation table <segregation_table>`
    :param float slice_thickness: Thickness of the NP.
    :param float nuclear_radius: Average nuclear radius.
    :param float coverage_q: The quantile at which to calculate window coverage. \
          The default value of 0.2 returns the coverage exceeded by 80% of \
          genomic windows.
    :param list skip_chroms: List of chromosome names to exclude.
    :param float genome_size: Total size of the genome (if None, calculate from \
          the segregation table.
    :param int plexity: Number of nuclear profiles in each sample.
    :param bool only_visible: If true, don't include genomic windows that are \
          never detected across the dataset in calculation. This is the default \
          behaviour.
    :returns: False if either of the checks fail, otherwise True.
    """

    if skip_chroms is None:
        skip_chroms = []

    pass_qc = True

    window_coverage = get_segregation_coverage_qc(segregation_table, skip_chroms,
                                only_visible, coverage_q)

    if window_coverage >= 20:
        print('Coverage check passed. 80% of genomic windows are detected at '
              'least {:.0f} times'.format(window_coverage))
    else:
        print('Coverage check failed. 80% of genomic windows are detected at '
              'least {:.0f} times. This value should be at least 20, try increasing '
              'the genomic window size for this dataset.'.format(window_coverage))
        pass_qc = False

    detection_efficiency = get_detection_efficiency(
        segregation_table, slice_thickness, nuclear_radius,
        skip_chroms, genome_size, plexity, only_visible)

    if detection_efficiency > 0.7:
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
              '{:.1%}. This value should be at least 70%.'.format(
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
