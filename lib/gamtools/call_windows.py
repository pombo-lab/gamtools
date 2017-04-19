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

from __future__ import print_function
import sys
import os
from bisect import bisect_left

import pandas as pd
import numpy as np
from scipy.stats import scoreatpercentile, nbinom, norm
from scipy.optimize import fmin

# define plt as a global so that we can import
# matplotlib later if we need it
plt = None


def cumulative_neg_binom(x, n, p):
    """
    Get the cumulative probability distribution for a negative
    binomial, over an array of points x that are log-scaled.

    :param x: Points at which to calculate probability density
    :type x: :class:`~numpy.ndarray`
    :param int n: Number of trials (see :data:`~scipy.stats.nbinom`)
    :param int p: Probability of success (see :data:`~scipy.stats.nbinom`)
    :returns: Cumulative probability distribution over x
    """

    # x is in log, so transform it back
    x = list(10 ** x[1:])

    # Add the point 0.0
    x = [0.0] + x

    return nbinom.cdf(x, n, p)


def un_cumulative(x):
    """
    Take a cumulative probability distribution x and extract
    the underlying probability density function.

    :param x: A cumulative probability distribution.
    :type x: :class:`~numpy.ndarray`
    :returns: Probability density distribution
    """

    last = None
    new_x = []
    for i in x:
        if last is None:
            last = i
            continue
        new_x.append(i - last)
        last = i
    return np.array(new_x)


def sum_to_1(x):
    """
    Normalize an array x such that the sum of all values is 1.

    :param x: Input array
    :type x: :class:`~numpy.ndarray`
    :returns: Normalized version of x
    """

    return x / sum(x)


def neg_binomial(x, n, p):
    """
    Get the probability density for a negative binomial distribution.

    :param x: Edges of bins within which to calculate probability density
    :type x: :class:`~numpy.ndarray`
    :param int n: Number of trials (see :data:`~scipy.stats.nbinom`)
    :param int p: Probability of success (see :data:`~scipy.stats.nbinom`)
    :returns: Probability density over x. Return \
            array has one less value than x, as the first value of \
            the return array is the probability density between x[0] \
            and x[1].
    """

    bin_y = cumulative_neg_binom(x, n, p)
    bin_y = un_cumulative(bin_y)

    return sum_to_1(bin_y)


def normal(x, loc, scale):
    """
    Get the probability density for a normal distribution.

    :param x: Edges of bins within which to calculate probability density
    :type x: :class:`~numpy.ndarray`
    :param int loc: Mean of the distribution (see :data:`~scipy.stats.norm`)
    :param int scale: Standard deviation of the distribution \
            (see :data:`~scipy.stats.norm`)
    :returns: Probability density over x. Return \
            array has one less value than x, as the first value of \
            the return array is the probability density between x[0] \
            and x[1].
    """

    norm_y = norm.cdf(x, loc, scale)
    norm_y = un_cumulative(norm_y)
    return sum_to_1(norm_y)


def n_binom_plus_log_normal(params, x):
    """
    Composite probability density function for the sum of a lognormal
    distribution and a negative binomial distribution.

    params is a tuple of all parameters for both underlying distributions:

    params = bin_n, bin_p, nm_delta, nm_scale, size

    Parameters for the negative binomial distribution are:

    bin_n = Negative binomial number of trials (see :data:`~scipy.stats.nbinom`)
    bin_p = Negative binomial probability of success (see :data:`~scipy.stats.nbinom`)

    Parameters for the lognormal distribution are:

    nm_delta = Absolute difference between the mean of the lognormal distribution
    and the mean of the negative binomial distribution
    nm_scale = Standard deviation of the lognormal distribution

    The lognormal is parameterized in this particular way because
    we don't want any solutions where the mean of the lognormal is
    less than the mean of the negative binomial. By parameterizing the function
    such that the position of the lognormal is given as a distance from the
    mean of the negative binomial (nm_delta), we can impose that nm_delta
    is always treated as positive.

    The final parameter, size, gives the ratio between the two underlying
    probability distributions.

    :param x: Edges of bins within which to calculate probability density
    :type x: :class:`~numpy.ndarray`
    :param tuple params: Parameters of the composite function
    :returns: Cumulative probability distribution over x. Return \
            array has one less value than x, as the first value of \
            the return array is the probability density between x[0] \
            and x[1].
    """

    bin_n, bin_p, nm_delta, nm_scale, size = params

    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)

    nm_loc = np.log10(bin_mean) + np.abs(nm_delta)

    bin_y = neg_binomial(x, bin_n, bin_p)

    norm_y = normal(x, nm_loc, nm_scale)

    sum_y = (bin_y * (1. - abs(size))) + (norm_y * abs(size))

    return sum_y / sum(sum_y)


def get_fdr_threshold(x, fdr, threshold):
    """
    Threshold a false-discovery rate distribution

    :param x: Positions at which false discovery rate is given
    :type x: :class:`~numpy.ndarray`
    :param fdr: False discovery rate (i.e. false positives divided by \
            total positives) at each position x.
    :type fdr: :class:`~numpy.ndarray`
    :param int threshold: Maximum false discovery rate at which to threshold
    :returns: Minimum values of x where fdr falls below the specified threshold.
    """

    threshold_index = bisect_left(fdr[::-1], threshold)
    return x[::-1][threshold_index - 1]


def squared_difference(params, my_function, x, y):
    """
    Get the squared difference between the output of an arbitrary function
    called with arbitrary parameters, and the expected output.

    This function is useful for fitting probability distributions (e.g
    using scipy's :data:`~scipy.optimize.fmin` function.

    :param tuple params: Tuple of parameters to be passed to my_function
    :param func my_function: Arbitrary function. Must accept exactly two \
            parameters, a tuple of parameters params and a numpy array x \
            of positions, and must itself return a numpy array.
    :param x: Positions to calculated the value of the function
    :type x: :class:`~numpy.ndarray`
    :param y: Expected value of the function at each position x.
    :type y: :class:`~numpy.ndarray`
    """

    z = my_function(params, x)
    diff = [d**2 for d in z - y]
    return sum(diff)


def mask_x_by_z(x, z):
    """
    Given two arrays x and z of equal length, return a list of the
    values from x excluding any where the corresponding value of z is
    0.

    :param x: Output values are taken from this array
    :type x: :class:`~numpy.ndarray`
    :param z: x[i] is excluded from the output if z[i] == 0
    :type z: :class:`~numpy.ndarray`
    """

    assert len(x) == len(z)

    return [xi for xi, zi in zip(x, z) if zi != 0]


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


def filter_data(x, percentile, no_zeros=True):
    """Remove data from an array which is below a certain
    percentile value. Optionally, if no_zeros is specified,
    also remove any zeros from the array.

    If removing values would result in returning an empty array,
    do nothing.

    :param x: Output values are taken from this array
    :type x: :class:`~numpy.ndarray`
    :param float percentile: Percentile at which to remove values \
            (e.g. if percentile=95.0, only the top 5% of values \
            are retained).
    :param bool no_zeros: If True, also discard any values equal \
            to zero from the output array.
    :returns: New array contining values from x that pass the filter.
    """

    percentile_score = scoreatpercentile(x, percentile)
    less_than_percentile = list(x < percentile_score)

    if no_zeros:
        not_a_zero = x > 0

        # only keep points which are both less than percentile AND not a zero
        points_to_keep = list(map(all, list(zip(less_than_percentile, not_a_zero))))

    else:
        points_to_keep = less_than_percentile

    out_data = x[points_to_keep]

    if out_data.size:

        return out_data

    if no_zeros:

        return x[not_a_zero]

    return x


def threshold_n_binom(params, p_value, thresh_range=None):
    """
    Determine a p-value threshold for a composite negative binomial
    and lognormal distribution based only on the value of the negative
    binomial.

    :param tuple params: Tuple of parameters for a combined \
            negative/binomial (see :func:`~n_binom_plus_log_normal`)
    :param float p_value: P-value cut-off
    :param list thresh_range: Possible values to consider as a cut \
            off (default is 0-500).
    :returns: Position above which the integral of the negative \
            binomial is equal to the P-value cut-off.
    """

    if thresh_range is None:
        thresh_range = list(range(500))

    bin_n, bin_p, nm_delta, nm_scale, size = params

    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)

    cumulative_dist = nbinom.cdf(thresh_range, bin_n, bin_p)

    prob_dist = sum_to_1(un_cumulative(cumulative_dist))
    index = bisect_left(prob_dist[::-1], p_value)
    return thresh_range[::-1][index]


def fit_histogram(breaks, counts, initial_params=None):
    """
    Fit a composite negative binomial and lognormal distribution
    to a histogram of read coverage per genomic window by
    least-squares minimization using the nelder-mead simplex
    algorithm.

    :param breaks: Array containing histogram bin edges.
    :type breaks: :class:`~numpy.ndarray`
    :param counts: Array containing the number of windows \
            falling into each histogram bin.
    :type counts: :class:`~numpy.ndarray`
    :param tuple initial_params: Initial guess for parameters \
            of the composite distribution (see \
            :func:`~n_binom_plus_log_normal`).
    :returns: Tuple of parameters representing the best-fit \
            distribution.
    """

    if initial_params is None:
        initial_params = (6.89799811e-01, 5.08503431e-01,
                          2.69945316, 0.23982432,
                          0.15)

    old_stdout = sys.stdout
    with open(os.devnull, "w") as sys.stdout:
        params = fmin(squared_difference,
                      initial_params,
                      (n_binom_plus_log_normal, breaks,
                       counts / float(sum(counts))),
                      xtol=0.00001, ftol=0.00001)
    sys.stdout = old_stdout
    return params


def get_fit_x(breaks, counts):
    """
    Get the x-values for plotting a fit based on the
    composite distribution.

    :param breaks: Array containing histogram bin edges.
    :type breaks: :class:`~numpy.ndarray`
    :param counts: Array containing the number of windows \
            falling into each histogram bin.
    :type counts: :class:`~numpy.ndarray`
    :returns: Array of x-values for plotting.
    """

    step = (breaks[1] - breaks[0]) / 2.
    fit_x = mask_x_by_z(breaks[:-1] + step, counts)
    return fit_x


def plot_combined_signal_noise(breaks, counts, params):
    """
    Plot a composite negative binomial and lognormal
    distribution.

    :param breaks: Array containing histogram bin edges.
    :type breaks: :class:`~numpy.ndarray`
    :param counts: Array containing the number of windows \
            falling into each histogram bin.
    :type counts: :class:`~numpy.ndarray`
    :param tuple params: Parameters of the composite \
            distribution (see  :func:`~n_binom_plus_log_normal`).
    """
    global plt

    if plt is None:
        try:
            from matplotlib import pyplot
            plt = pyplot
        except ImportError:
            raise ImportError('Plotting requires matplotlib to be installed.')

    fit = sum(counts) * n_binom_plus_log_normal(params,
                                                breaks)

    fit_y = mask_x_by_z(fit, counts)

    fit_x = get_fit_x(breaks, counts)

    return plt.plot(fit_x, fit_y, 'ro')


def plot_lognormal(breaks, counts, params):
    """
    Given the parameters of a composite negative binomal
    and lognormal distribution, plot only the contribution
    of the lognormal.

    :param breaks: Array containing histogram bin edges.
    :type breaks: :class:`~numpy.ndarray`
    :param counts: Array containing the number of windows \
            falling into each histogram bin.
    :type counts: :class:`~numpy.ndarray`
    :param tuple params: Parameters of the composite \
            distribution (see  :func:`~n_binom_plus_log_normal`).
    """

    bin_n, bin_p, nm_delta, nm_scale, size = params

    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)

    nm_loc = np.log10(bin_mean) + np.abs(nm_delta)

    gauss_y = normal(breaks, nm_loc, nm_scale)
    gauss_y = gauss_y * sum(counts) * abs(size)
    gauss_y = mask_x_by_z(gauss_y, counts)

    fit_x = get_fit_x(breaks, counts)

    return plt.plot(fit_x, gauss_y, color='yellow')


def plot_binom(breaks, counts, params):
    """
    Given the parameters of a composite negative binomal
    and lognormal distribution, plot only the contribution
    of the negative binomial.

    :param breaks: Array containing histogram bin edges.
    :type breaks: :class:`~numpy.ndarray`
    :param counts: Array containing the number of windows \
            falling into each histogram bin.
    :type counts: :class:`~numpy.ndarray`
    :param tuple params: Parameters of the composite \
            distribution (see  :func:`~n_binom_plus_log_normal`).
    """

    bin_n, bin_p, nm_delta, nm_scale, size = params

    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)

    binom_y = neg_binomial(breaks, bin_n, bin_p)
    binom_y = binom_y * sum(counts) * (1. - abs(size))
    binom_y = mask_x_by_z(binom_y, counts)

    fit_x = get_fit_x(breaks, counts)

    return plt.plot(fit_x, binom_y, color='green')


def plot_legend(hist_patches,
                fit_patches,
                normal_patches,
                binom_patches,
                thresh_patch):
    """
    Make a legend for a plot of the composite fitting.
    """

    patches = (hist_patches[0],
               fit_patches[0],
               normal_patches[0],
               binom_patches[0],
               thresh_patch)

    labels = ('Coverage histogram',
              'Fitted values',
              'Signal distribution',
              'Noise distribution',
              'Threshold')

    plt.legend(patches, labels)


def prettify_plot(sample_name, fig):
    """
    Add titles and axis labels for a plot of the composite fitting.
    """

    fig.suptitle('Combined fit for {0}'.format(sample_name), fontsize=18)

    plt.ylabel('No of windows')
    plt.xlabel('No of reads')
    locs, labels = plt.xticks()
    labels = list(map(int, 10 ** locs))
    plt.xticks(locs, labels)


def plot_signal_and_noise_fitting(sample_name, fitting_results):
    """
    Make a plot of a composite negative binomial and lognormal
    distribution fit to a read coverage histogram.

    :param str sample_name: Name of the sample to plot.
    :param dict fitting_results: Dictionary of results returned \
            by :func:`~signal_and_noise_fitting`.
    """

    global plt

    if plt is None:
        try:
            from matplotlib import pyplot
            plt = pyplot
        except ImportError:
            raise ImportError('Plotting requires matplotlib to be installed.')

    fig = plt.figure(figsize=(16, 9))

    bar_width = (fitting_results['breaks'].max() -
                 fitting_results['breaks'].min()) / 50.
    hist_patches = plt.bar(
        fitting_results['breaks'][
            :-1],
        fitting_results['counts'],
        width=bar_width)

    fit_patches = plot_combined_signal_noise(fitting_results['breaks'],
                                             fitting_results['counts'],
                                             fitting_results['params'])

    normal_patches = plot_lognormal(fitting_results['breaks'],
                                    fitting_results['counts'],
                                    fitting_results['params'])

    binom_patches = plot_binom(fitting_results['breaks'],
                               fitting_results['counts'],
                               fitting_results['params'])

    thresh_patch = plt.axvline(
        np.log10(
            fitting_results['read_threshold']),
        color='black')

    plot_legend(hist_patches,
                fit_patches,
                normal_patches,
                binom_patches,
                thresh_patch)

    prettify_plot(sample_name, fig)


def fixed_threshold_fitting_func(read_threshold):
    """
    Factory function for creating fitting functions that return a pre-specified
    threshold for every sample.

    :param int read_threshold: Read threshold to use for every sample.
    :returns: Fitting function.
    """

    read_threshold = int(read_threshold)

    def fixed_threshold_fitting(sample_coverage_data):
        """
        Return a pre-specified read coverage threshold.

        :param sample_coverage_data: Array of read counts per \
                genomic window.
        :type sample_coverage_data: :class:`~numpy.ndarray`
        :returns: Dictionary of fitting results.
        """

        counts, breaks = np.histogram(
            np.log10(
                filter_data(
                    sample_coverage_data, 99.99)), bins=50)

        params = None

        return {'read_threshold': read_threshold,
                'counts': counts,
                'breaks': breaks,
                'params': params}

    return fixed_threshold_fitting


def signal_and_noise_fitting(sample_coverage_data):
    """
    Fit a composite negative binomial and lognormal
    distribution to set of read counts per genomic window.

    :param sample_coverage_data: Array of read counts per \
            genomic window.
    :type sample_coverage_data: :class:`~numpy.ndarray`
    :returns: Dictionary of fitting results.
    """

    counts, breaks = np.histogram(
        np.log10(
            filter_data(
                sample_coverage_data, 99.99)), bins=50)

    params = fit_histogram(breaks, counts)

    read_threshold = threshold_n_binom(params, 0.001)

    return {'read_threshold': read_threshold,
            'counts': counts,
            'breaks': breaks,
            'params': params}


def plot_fitting_and_save(sample_name, fitting_folder, sample_fitting_data):
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

    plot_signal_and_noise_fitting(sample_name, sample_fitting_data)

    plot_path = os.path.join(fitting_folder,
                             '{0}_fit.png'.format(safe_sample_name))

    if not os.path.isdir(fitting_folder):
        os.makedirs(fitting_folder)

    plt.savefig(plot_path)

    plt.close()


def do_coverage_thresholding(coverage_data, fitting_folder, fitting_function):
    """
    Given a :ref:`read coverage table <read_coverage_table>`, fit a composite negative binomial and
    lognormal distribution to the read coverage histogram for each individual
    NP. Using this distribution, determine a read coverage threshold to
    distinguish positive (signal) windows from negative (noise) windows.
    Finally, use the determined thresholds to convert the
    :ref:`read coverage table <read_coverage_table>` to a
    :ref:`segregation table <segregation_table>` giving the
    locations of positive windows for each NP.

    :param coverage_data: :ref:`read coverage table <read_coverage_table>`
    :param str fitting_folder: Path to the folder in which to save \
            images of each individual fit.
    :param func fitting_function: Function to use to determine a read \
            threshold for each NP.
    :returns: segregation_table (a :ref:`segregation table <segregation_table>`) and \
            fitting_data (a :class:`~pandas.DataFrame` of fitting \
            parameters for each NP).
    """

    fitting_list = []

    segregation_table = coverage_data.copy()

    for i, sample_name in enumerate(coverage_data.columns):

        sample_coverage_data = coverage_data.loc[:, sample_name]

        sample_fitting_data = fitting_function(sample_coverage_data)

        if fitting_folder is not None:

            plot_fitting_and_save(
                sample_name,
                fitting_folder,
                sample_fitting_data)

        above_threshold = coverage_data.loc[
            :, sample_name] > sample_fitting_data['read_threshold']

        segregation_table.loc[:, sample_name] = above_threshold.astype(int)

        fitting_list.append([sample_name] + list(sample_fitting_data.values()))

        print(
            '{0} done (number {1}) threshold {2}'.format(
                sample_name,
                i + 1,
                sample_fitting_data['read_threshold']),
            file=sys.stderr)

    fitting_data = pd.DataFrame(
        fitting_list, columns=['Sample'] + list(sample_fitting_data.keys())
    )

    return segregation_table, fitting_data


def threshold_file(input_file, output_file,
                   fitting_folder, fitting_details_file, fitting_function):
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
    :param func fitting_function: Function to use for thresholding each NP.
    """

    coverage_data = pd.read_csv(input_file,
                                delim_whitespace=True, index_col=[0, 1, 2])

    segregation_matrix, fitting_data = do_coverage_thresholding(
        coverage_data, fitting_folder, fitting_function)

    segregation_matrix.to_csv(output_file, index=True, sep='\t')

    if fitting_details_file is not None:
        fitting_data.to_csv(fitting_details_file, sep='\t', index=False)


def threshold_from_args(args):
    """
    Wrapper function to pass arguments to threshold_file from argparse.
    """

    if args.macs:
        raise NotImplementedError('Thresholding using macs is not supported')

    else:
        threshold_file(args.coverage_file, args.output_file,
                       args.fitting_folder, args.details_file,
                       args.fitting_function)
