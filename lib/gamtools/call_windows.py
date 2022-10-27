"""
=======================
The call_windows module
=======================


Calling positive windows
========================

The call windows module contains code for identifying positive
window from GAM sequencing data.

.. _read_coverage_table:

Read Coverage Tables
====================

Read coverage tables are tables that give the number of sequencing reads from
each NP that map to each genomic region. Each column gives the read coverage
for a single NP, whilst each row gives the number of reads mapping to a
specific genomic region. In order to proceed further with analyzing the output
of a GAM experiment, read coverage tables must be converted to
:ref:`segregation_table` by identifying windows which were truly present in
each NP and discarding regions where a low number of mapped reads are more
likely to stem from background sequencing noise than from genomic DNA that was
truly present in the original NP.

"""

# pylint: disable=invalid-name,unused-variable,global-statement

import sys
import os

import pandas as pd
import numpy as np
from scipy.stats import mode
from . import segregation

# Don't throw an error if matplotlib is not installed unless we try to use it
try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = None

def erode(coverage_data, prob):
    """
    Probabilistically erode a
    :ref:`read coverage table <read_coverage_table>`, to generate a new table
    with less total read coverage. Each entry in the table, representing a
    number of reads mapping to a given window, is replaced by a random value
    drawn from a binomial distribution where each of the original reads has a
    certain probability of detection.

    :param coverage_data: Read coverage table to erode.
    :param float prob: Amount by which to erode the original table (e.g. if \
            the original table had 1000 reads and prob is 0.7, the new table \
            will have approximately 700 reads.
    :returns: New :ref:`read coverage table <read_coverage_table>` \
            with eroded values
    """

    return coverage_data.astype(np.int32).apply(
        lambda x: np.random.binomial(x, prob))


def fixed_threshold_fitting_func(coverage_threshold):
    """
    Factory function for creating fitting functions that return a pre-specified
    threshold for every sample.

    :param int coverage_threshold: Read threshold to use for every sample.
    :returns: Fitting function.
    """

    coverage_threshold = float(coverage_threshold)

    def fixed_threshold_fitting(sample_coverage_data): # pylint: disable=unused-argument
        """
        Return a pre-specified read coverage threshold.

        :param sample_coverage_data: Array of read counts per \
                genomic window.
        :type sample_coverage_data: :class:`~numpy.ndarray`
        :returns: Dictionary of fitting results.
        """

        params = None

        return {'coverage_threshold': coverage_threshold,
                'params': params}

    return fixed_threshold_fitting



def pct_orphan_windows(positive_windows):
    """Given an array of windows, calculate the percentage of orphan
    windows (positive windows not neighbouring other positive windows).

    :param positive_windows: Array representing positive or negative windows.
    :returns: Percentage of positive windows not bordering other positive \
    windows.
    """

    # If there are no positive windows we cannot calculate a percentage
    # of orphan windows
    if positive_windows.sum() == 0:
        return np.Inf

    # Orphan windows are positive with two negative neighbours (0, 1, 0)
    window_filter = np.array([[0], [1], [0]])

    # construct Rolling 3-bin window over input data
    rolling_window = np.vstack(
        (positive_windows[:-2], positive_windows[1:-1], positive_windows[2:]))

    # check which positions in the rolling window match the filter
    orphan_windows = np.all(rolling_window == window_filter, axis=0)

    total_orphans = orphan_windows.sum()

    # Check whether the first and last windows are orphans
    if (positive_windows[0] == 1) and (positive_windows[1] == 0):
        total_orphans += 1

    if (positive_windows[-1] == 1) and (positive_windows[-2] == 0):
        total_orphans += 1

    return float(total_orphans) / positive_windows.sum()


def infer_single_read_coverage(coverage_table):
    """
    Infer the coverage of a single read given the read coverage table.
    Specifically, find the most frequent non-zero entry in the coverage matrix and
    check that the value is not exactly double any of the other 20 most frequent
    entries (in case the most frequent value does not represent one read). Return
    the inferred coverage from one read multiplied by the inferred window size.

    :param sample_coverage_data: Array of bp coverage per \
        genomic window.
    :returns: Inferred coverage of one single read.
    """

    cov_values, cov_counts = np.unique(coverage_table, return_counts=True)
    top_20_vals = cov_values[cov_counts.argsort()[-20:]]
    top_20_vals = top_20_vals[top_20_vals > 0]
    top_val = top_20_vals[-1]

    half_top_val = np.logical_and(
        (top_20_vals / top_val) > 0.45,
        (top_20_vals / top_val) < 0.55)

    if np.any(half_top_val):
        top_val = top_20_vals[half_top_val]

    return top_val


def orphan_windows_fitting_func(percentile_str, clip=True):
    """
    Factory function for creating fitting functions that find the threshold
    giving the minimum percentage of orphan windows.

    :param str percentile_str: Upper and lower bounds for the percentiles of \
      the coverage distribution to search for a threshold, for example \
      '10,80' would force a threshold that is higher than 10% of the windows \
      with non-zero coverage, but lower than 80%.
    :param bool clip: Whether to clip the threshold to a minimum of one read \
      (i.e. make sure windows need more than one mapping read to be \
      considered positive)
    :returns: Fitting function.
    """

    bottom_percentile, top_percentile = (int(pct) for pct in percentile_str.split(','))

    def orphan_windows_fitting(sample_coverage_data):
        """ Find a coverage threshold that gives the minimum percentage of orphan
        windows (positive windows not neighbouring other positive windows).

        :param sample_coverage_data: Array of bp coverage per \
                genomic window.
        :type sample_coverage_data: :class:`~numpy.ndarray`
        :returns: Dictionary of fitting results.
        """

        nonzero_coverage = sample_coverage_data[sample_coverage_data > 0]

        thresholds = [np.percentile(nonzero_coverage, pct)
                      for pct in range(bottom_percentile, top_percentile + 1)]
        orphan_windows = []

        for threshold in thresholds:

            positive_windows = sample_coverage_data > threshold

            orphan_windows.append(
                pct_orphan_windows(positive_windows))

        optimal_threshold = thresholds[np.argmin(orphan_windows)]

        if clip:

            single_read_estimate = mode(nonzero_coverage)[0][0]

            optimal_threshold = np.clip(optimal_threshold, single_read_estimate, None)

        return {'coverage_threshold': optimal_threshold}

    return orphan_windows_fitting


def prettify_plot(sample_name, fig):
    """
    Add titles and axis labels for a plot of the composite fitting.
    """

    fig.suptitle('Combined fit for {0}'.format(sample_name), fontsize=18)

    plt.ylabel('No of windows')
    plt.xlabel('Average depth')
    locs, labels = plt.xticks()
    labels = ['{:.4f}'.format(10 ** loc) for loc in locs]
    plt.xticks(locs, labels)
    plt.legend()


def plot_coverage_histogram(sample_name, sample_coverage_data, sample_fitting_data):
    """
    Plot the coverage histogram and threshold for a given sample.

    :param str sample_name: Name of the sample to use as plot title
    :param array sample_coverage_data: Array of windows with read coverages
    :param dict sample_fitting_data: Dictionary containing a \
      'coverage_threshold' key
    """

    if plt is None:
        raise ImportError('Plotting requires matplotlib to be installed.')

    fig = plt.figure(figsize=(16, 9))

    counts, breaks = np.histogram(
        np.log10(sample_coverage_data[sample_coverage_data > 0]), bins=50)

    bar_width = (breaks.max() - breaks.min()) / 50

    plt.bar(breaks[:-1], counts, width=bar_width)
    plt.axvline(np.log10(sample_fitting_data['coverage_threshold']) + (bar_width / 2),
                color='red', label='Threshold')

    prettify_plot(sample_name, fig)


def plot_fitting_and_save(sample_name, fitting_folder, sample_coverage_data, sample_fitting_data):
    """
    Make a plot of a composite negative binomial and lognormal
    distribution fit to a read coverage histogram, then save the
    plot to a file.

    :param str sample_name: Name of the sample to plot.
    :param str fitting_folder: Path to the folder in which to save \
            the image.
    :param dict fitting_results: Dictionary of results returned \
            by :func:`~signal_and_noise_fitting`.
    """

    safe_sample_name = os.path.basename(sample_name)

    plot_coverage_histogram(sample_name, sample_coverage_data, sample_fitting_data)

    plot_path = os.path.join(fitting_folder,
                             '{0}_fit.png'.format(safe_sample_name))

    if not os.path.isdir(fitting_folder):
        os.makedirs(fitting_folder)

    plt.savefig(plot_path)

    plt.close()


def do_coverage_thresholding(coverage_data, fitting_folder, fitting_function,
			     min_reads=None):
    """
    Given a :ref:`read coverage table <read_coverage_table>`, extract the read
    coverage histogram for each individual NP. Pass this distribution to a
    function to determine a read coverage threshold to distinguish positive
    (signal) windows from negative (noise) windows.  Finally, use the
    determined thresholds to convert the :ref:`read coverage table
    <read_coverage_table>` to a :ref:`segregation table <segregation_table>`
    giving the locations of positive windows for each NP.

    :param coverage_data: :ref:`read coverage table <read_coverage_table>`
    :param str fitting_folder: Path to the folder in which to save \
            images of each individual fit.
    :param func fitting_function: Function to use to determine a read \
            threshold for each NP.
    :param min_reads: Minimum number of reads required to call a window as \
            positive.
    :type min_reads: int or None
    :returns: segregation_table (a :ref:`segregation table <segregation_table>`) and \
            fitting_data (a :class:`~pandas.DataFrame` of fitting \
            parameters for each NP).
    """

    fitting_list = []

    segregation_table = coverage_data.copy()

    if min_reads is not None:
        inferred_single_read_cov = infer_single_read_coverage(coverage_data)
        min_cov = min_reads * inferred_single_read_cov
    else:
        min_cov = 0

    for i, sample_name in enumerate(coverage_data.columns):

        sample_coverage_data = coverage_data.loc[:, sample_name]

        sample_fitting_data = fitting_function(sample_coverage_data)

        sample_fitting_data['coverage_threshold'] = max(
            sample_fitting_data['coverage_threshold'], min_cov)

        if fitting_folder is not None:

            plot_fitting_and_save(
                sample_name,
                fitting_folder,
                sample_coverage_data,
                sample_fitting_data)

        above_threshold = coverage_data.loc[
            :, sample_name] > sample_fitting_data['coverage_threshold']

        segregation_table.loc[:, sample_name] = above_threshold.astype(int)

        fitting_list.append([sample_name] + list(sample_fitting_data.values()))

        print(
            '{0} done (number {1}) threshold {2}'.format(
                sample_name,
                i + 1,
                sample_fitting_data['coverage_threshold']),
            file=sys.stderr)

    fitting_data = pd.DataFrame(
        fitting_list, columns=['Sample'] + list(sample_fitting_data.keys())
    )

    return segregation_table, fitting_data


def threshold_file(input_file, output_file, #pylint: disable=too-many-arguments
                   fitting_folder, fitting_details_file,
                   fitting_function, min_reads=None):
    """
    Read a :ref:`read coverage table <read_coverage_table>` file and threshold
    each NP, to generate a :ref:`segregation table <segregation_table>`.

    :param str input_file: Path to \
            :ref:`read coverage table <read_coverage_table>`.
    :param str output_file: Path to save \
            :ref:`segregation table <segregation_table>`.
    :param fitting_folder: Path to save plots of each individual fit \
            (or None to skip saving fitting images).
    :type fitting_folder: str or None
    :param fitting_details_file: Path to save table of parameters for each \
            fit (or None to skip saving fitting table).
    :type fitting_details_file: str or None
    :param min_reads: Minimum number of reads required to call a windows \
            as positive.
    :type min_reads: int or None
    :param func fitting_function: Function to use for thresholding each NP.
    """

    coverage_data = pd.read_csv(input_file,
                                delim_whitespace=True, index_col=[0, 1, 2])

    segregation_matrix, fitting_data = do_coverage_thresholding(
        coverage_data, fitting_folder, fitting_function, min_reads)

    segregation_matrix.to_csv(output_file, index=True, sep='\t')

    if fitting_details_file is not None:
        fitting_data.to_csv(fitting_details_file, sep='\t', index=False)


def merge_reads_coverage(targets, dependencies):
    """Merge multiple files containing coverage per window for single libraries
    into one big table containing one column for each input file.

    :param targets: An iterable containing one output path for saving the \
            merged table
    :param dependencies: An iterable containing the paths for all of the \
            input files
    """

    (output_file, ) = targets
    input_coverage_paths = list(dependencies)

    first_coverage_path = input_coverage_paths.pop()
    merged_coverage_table = segregation.open_segregation(first_coverage_path)

    for next_coverage_path in input_coverage_paths:
        merged_coverage_table = pd.merge(
            merged_coverage_table,
            segregation.open_segregation(next_coverage_path),
            left_index=True, right_index=True)

    merged_coverage_table.to_csv(output_file, sep='\t')


def threshold_from_args(args):
    """
    Wrapper function to pass arguments to threshold_file from argparse.
    """

    threshold_file(args.coverage_file, args.output_file,
                   args.fitting_folder, args.details_file,
                   args.fitting_function, args.min_reads)
