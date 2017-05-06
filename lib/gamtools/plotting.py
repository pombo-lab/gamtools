"""
===================
The plotting module
===================

The plotting module contains functions for plotting the read coverage
and the positive windows identified in a given NP.

"""

#pylint: disable=invalid-name,too-many-arguments,too-many-locals

import numpy as np
import pandas as pd
from .utils import DelayedImportError

failed_import = DelayedImportError(message=None)
failed_packages = []

try:
    from pybedtools import Interval
except ImportError:
    Interval = failed_import
    failed_packages.append('pybedtools')

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = failed_import
    failed_packages.append('matplotlib')

try:
    from metaseq import genomic_signal
except ImportError:
    genomic_signal = failed_import
    failed_packages.append('metaseq')

if failed_packages:

    if len(failed_packages) > 1:
        package_list = ', '.join(failed_packages[:-1])
        package_list += ' and '
        package_list += failed_packages[-1]
    else:
        package_list = failed_packages[0]

    pip_call = ' '.join(failed_packages)

    base_message = ('GAMtools plotting commands require {package_list} to'
                    ' be installed. Try to install them by running '
                    '"pip install {pip_call}"')

    failed_import.message = base_message.format(package_list=package_list,
                                                pip_call=pip_call)


def chunk_genomic_signal(
        genomic_signal_obj,
        interval,
        bins=500,
        chunk_size=20000000):
    """Retrieve intervals from a MetaSeq genomic_signal object in chunks, instead
    of all at once.

    :param genomic_signal_obj: Any MetaSeq genomic_signal object to query
    :param interval: Region of the genome to query
    :type interval: PyBedTools Interval object
    :param int bins: Total number of bins to return
    :param int chunk_size: Size in bp of each chunk

    :returns: array of genomic starts for each bin, \
            array of BigBed signal values for each bin
    """

    def _internal_chunker():
        for start in range(interval.start, interval.end, chunk_size):

            end = start + chunk_size
            if end > interval.end:
                end = interval.end

            length = end - start
            chunk_bins = int((float(length) / interval.length) * bins)

            new_interval = Interval(interval.chrom, start, end)

            try:
                chunk_x, chunk_y = genomic_signal_obj.local_coverage(
                    new_interval, bins=chunk_bins)
            except ValueError:
                chunk_x = np.linspace(start, end, chunk_bins)
                chunk_y = np.zeros_like(chunk_x)

            yield chunk_x, chunk_y

    x_chunks, y_chunks = list(zip(*list(_internal_chunker())))

    x, y = np.concatenate(x_chunks), np.concatenate(y_chunks)

    return x, y


def open_sizes_file(sizes_path):
    """Open a chromosome sizes file of the sort provided by UCSC

    :param str sizes_path: Path to chromosome sizes file

    :returns: Pandas DataFrame listing chromosome names and sizes
    """

    chrom_sizes = pd.read_csv(
        sizes_path,
        delim_whitespace=True,
        header=None,
        names=[
            'chrom',
            'size'],
        index_col=0)

    return chrom_sizes


def assign_chroms_to_rows(chrom_sizes):
    """Assign chromosomes to rows of a plot, such that each row covers 18%
    of the total genome.

    :param chrom_sizes: Table of chromosome names and sizes
    :type chrom_sizes: Pandas DataFrame

    :returns: List of lists of chromosome names by row, \
            list of lists of chromosome sizes by row
    """

    row_tot = 0.
    rows = []
    current_row = []
    chrom_pct_of_genome = list(
        chrom_sizes['size'] *
        100 /
        chrom_sizes['size'].sum())
    for i, c in enumerate(chrom_pct_of_genome):
        current_row.append(i + 1)
        row_tot += c

        # Once the current row covers more than 18% of the genome, start a new
        # one
        if row_tot > 18.:
            row_tot = 0.
            rows.append(current_row)
            current_row = []

    chrom_names_by_row = [['chr{}'.format(c) for c in row] for row in rows]
    chrom_sizes_by_row = [[int(chrom_sizes.loc[n]) for n in row]
                          for row in chrom_names_by_row]

    return chrom_names_by_row, chrom_sizes_by_row


def parse_sizes_file(sizes_path):
    """Open a chromosome sizes file and assign chromosomes to rows of a plot,
    such that each row covers 18% of the total genome.

    :param str sizes_path: Path to chromosome sizes file

    :returns: List of lists of chromosome names by row, \
            list of lists of chromosome sizes by row
    """

    chrom_sizes = open_sizes_file(sizes_path)

    chrom_names_by_row, chrom_sizes_by_row = assign_chroms_to_rows(chrom_sizes)

    return chrom_names_by_row, chrom_sizes_by_row


def get_row_pct(row_sizes):
    """Get the span of each chromosome as a percentage of the span of the largest row

    :param list row_sizes: List of lists of chromosome sizes by row

    :returns: List of lists of chromosome sizes as percentages by row
    """

    biggest_row = max([sum(row) for row in row_sizes])
    row_pcts = [[float(val) / biggest_row for val in row] for row in row_sizes]
    return row_pcts


def plot_chrom_sig(signal_obj, ax, interval, color, bins=500):
    """Plot the read coverage for a given interval.

    Automatically scales the y-axis based on the 98th percentile of signal.

    :param signal_obj: Any MetaSeq genomic_signal object
    :param ax: Matplotlib figure axis to plot to
    :param interval: PyBedtools Interval to plot
    :param color: Color to use for plotting (any format accepted by matplotlib).
    :param int bins: Number of bins to average over
    """

    arr = chunk_genomic_signal(signal_obj, interval, bins=bins)

    top = np.ceil(np.percentile(arr[1], 98.))
    if top < 2:
        top = 2

    ax.fill_between(arr[0], arr[1], color=color)
    ax.set_ylim(1, top)
    ax.set_xlim(arr[0][0], arr[0][-1])
    ax.axis('off')


def plot_chrom_seg(signal_obj, ax, interval, color, bins=500):
    """Plot the positive windows for a given interval.

    Automatically scales the y-axis between 0 and 1

    :param signal_obj: Any MetaSeq genomic_signal object
    :param ax: Matplotlib figure axis to plot to
    :param interval: PyBedtools Interval to plot
    :param color: Color to use for plotting (any format accepted by matplotlib).
    :param int bins: Number of bins to average over
    """

    arr = chunk_genomic_signal(signal_obj, interval, bins=bins)

    ax.fill_between(arr[0], arr[1], color=color)
    ax.set_ylim(0.2, 1.2)
    ax.set_xlim(arr[0][0], arr[0][-1])
    ax.axis('off')


def plot_both(objects, axes, name, end, color, bins=500):
    """Plot the read coverage and the positive window calls for a given chromosome
    of a single NP.

    :param objects: (read coverage object, positive window object) \
            - both MetaSeq genomic_signal objects
    :param axes: (read coverage axis, positive window axis) \
            - both Matplotlib figure axis
    :param name: Name of chromosome to plot
    :param end: Size of chromosome to plot
    :param color: Color to use for plotting (any format accepted by matplotlib).
    :param int bins: Number of bins to average over
    """

    i = Interval(name, 0, end)

    sig_obj, seg_obj = objects
    sig_ax, seg_ax = axes

    plot_chrom_sig(sig_obj, sig_ax, i, color, bins)
    plot_chrom_seg(seg_obj, seg_ax, i, color, bins)


def plot_genome(
        axis_sizes,
        row_names,
        row_sizes,
        read_cov_obj,
        pos_window_obj,
        out_file):
    """Plot the read coverage and the positive window calls of a single NP
    for every chromosome.

    :param axis_sizes: List of lists. Each row is a list of axis sizes, \
            one for each chromosome that should be plotted on the row
    :param row_names: List of lists. Each row is a list of chromosome names, \
            one for each chromosome that should be plotted on the row
    :param row_sizes: List of lists. Each row is a list of chromosome sizes, \
            one for each chromosome that should be plotted on the row
    :param read_cov_obj: MetaSeq genomic_signal object giving read coverage
    :param pos_window_obj: MetaSeq genomic_signal object giving positive windows
    :param out_file: Path to save image to.
    """

    fig = plt.figure(figsize=(28, 5))

    axes = []
    row_axes = []

    curr_pos = 0
    for rnum, row in enumerate(axis_sizes):
        for size in row:
            sig_axis = plt.subplot2grid(
                (20, 30), (rnum * 4, curr_pos), colspan=size, rowspan=3)
            seg_axis = plt.subplot2grid(
                (20, 30), ((rnum * 4) + 3, curr_pos), colspan=size)
            row_axes.append((sig_axis, seg_axis))
            curr_pos += size
        axes.append(row_axes)
        row_axes = []
        curr_pos = 0

        i = 0
    for p, row in enumerate(axes):
        for q, ax in enumerate(row):
            if i % 2:
                col = '#2ca25f'
            else:
                col = '#e34a33'

            plot_both((read_cov_obj, pos_window_obj), ax,
                      row_names[p][q], row_sizes[p][q], color=col, bins=10000)
            i += 1

    max_ylim = max([max([ax[0].get_ylim()[1] for ax in row]) for row in axes])

    for row in axes:
        for ax in row:
            ax[0].set_ylim(1, max_ylim)

    fig.savefig(out_file)


def plot_np(bigwig_file, bigbed_file, sizes_file, output_file):
    """Plot the read coverage and the positive window calls of a single NP
    for every chromosome.

    :param bigwig_file: Path to bigwig file for NP read coverage
    :param bigbed_file: Path to bigbed file for NP positive windows
    :param sizes_file: Path to chromosome sizes file
    :param output_file: Path to save image to
    """

    row_names, row_sizes = parse_sizes_file(sizes_file)
    row_pcts = get_row_pct(row_sizes)
    axis_sizes = [[int(round(val * 30)) for val in row] for row in row_pcts]

    read_coverage = genomic_signal(bigwig_file, 'bigwig')
    positive_windows = genomic_signal(bigbed_file, 'bed')

    plot_genome(
        axis_sizes,
        row_names,
        row_sizes,
        read_coverage,
        positive_windows,
        output_file)


def plot_np_from_args(args):
    """Wrapper function to call plot_np from argparse"""

    plot_np(args.bigwig_file, args.bed_file,
            args.genome_file, args.output_file)
