import numpy as np
import pandas as pd

try:
    from metaseq import genomic_signal
except ImportError:
    raise ImportError(
        "The 'plot_np' command requires metaseq to be installed "
        "(try installing it using 'pip install metaseq').")

try:
    from pybedtools import Interval
except ImportError:
    raise ImportError(
        "The 'plot_np' command requires pybedtools to be installed "
        "(try installing it using 'pip install pybedtools').")

from matplotlib import pyplot as plt

def chunk_bigbed(bigbed, interval, bins=500, chunk_size=20000000):
    def _internal_chunker():
        for start in range(interval.start, interval.end, chunk_size):

            end = start + chunk_size
            if end > interval.end:
                end = interval.end

            length = end - start
            chunk_bins = int((float(length) / interval.length) * bins)

            new_interval = Interval(interval.chrom, start, end)

            chunk_x, chunk_y = bigbed.local_coverage(new_interval, bins=chunk_bins)

            yield chunk_x, chunk_y

    x_chunks, y_chunks = zip(*list(_internal_chunker()))

    x, y = np.concatenate(x_chunks), np.concatenate(y_chunks)

    return x, y

def parse_sizes_file(sizes_path):
    chrom_sizes = pd.read_csv(sizes_path, delim_whitespace=True,
                              header=None, names=['chrom', 'size'], index_col=0)

    row_tot = 0.
    rows = []
    current_row = []
    for i, c in enumerate(list(chrom_sizes['size'] * 100 / chrom_sizes['size'].sum())):
        current_row.append(i+1)
        row_tot += c
        if row_tot > 18.:
            row_tot = 0.
            rows.append(current_row)
            current_row = []
    names = [ map(lambda c:'chr'+str(c), row) for row in rows ]
    sizes = [ map(lambda n:int(chrom_sizes.loc[n]), row) for row in names ]
    return names, sizes

def get_row_pct(row_sizes):
    biggest_row = max([ sum(row) for row in row_sizes ])
    row_pcts = [ [ float(val) / biggest_row for val in row ] for row in row_sizes ]
    return row_pcts

def plot_chrom_sig(signal_obj, ax, interval, color, bins=500):

    arr = chunk_bigbed(signal_obj, interval, bins=bins)

    top = np.ceil(np.percentile(arr[1], 98.))
    ax.fill_between(arr[0], arr[1], color=color)
    ax.set_ylim(1,top)
    ax.set_xlim(arr[0][0], arr[0][-1])
    ax.axis('off')

def plot_chrom_seg(seg_obj, ax, interval, color, bins=500):

    arr = chunk_bigbed(seg_obj, interval, bins=bins)

    top = np.ceil(np.percentile(arr[1], 98.))
    ax.fill_between(arr[0], arr[1], color=color)
    ax.set_ylim(1,top)
    ax.set_xlim(arr[0][0], arr[0][-1])
    ax.axis('off')

def plot_both(objects, axes, name, end, color, bins=500):
    i = Interval(name, 0, end)

    sig_obj, seg_obj = objects
    sig_ax, seg_ax = axes

    plot_chrom_sig(sig_obj, sig_ax, i, color, bins)
    plot_chrom_sig(seg_obj, seg_ax, i, color, bins)

def plot_genome(axis_sizes, row_names, row_sizes, g, bb, out_file):

    fig = plt.figure(figsize=(28,5))

    axes = []
    row_axes = []

    curr_pos = 0
    for rnum, row in enumerate(axis_sizes):
        for size in row:
            sig_axis = plt.subplot2grid((20,30), (rnum*4, curr_pos), colspan=size, rowspan=3)
            seg_axis = plt.subplot2grid((20,30), ((rnum*4)+3, curr_pos), colspan=size)
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

            plot_both((g, bb), ax, row_names[p][q], row_sizes[p][q], color=col)
            i += 1

    max_ylim = max([ max([ ax[0].get_ylim()[1] for ax in row ]) for row in axes ])

    for row in axes:
        for ax in row:
            ax[0].set_ylim(1, max_ylim)

    fig.savefig(out_file)


def plot_np(bigwig_file, bigbed_file, sizes_file, output_file):

    row_names, row_sizes = parse_sizes_file(sizes_file)
    row_pcts = get_row_pct(row_sizes)
    axis_sizes = [ [ int(round(val * 30)) for val in row ] for row in row_pcts ]

    read_coverage = genomic_signal(bigwig_file, 'bigwig')
    positive_windows = genomic_signal(bigbed_file, 'bigbed')

    plot_genome(axis_sizes, row_names, row_sizes, read_coverage, positive_windows, output_file)


def plot_np_from_args(args):

    plot_np(args.bigwig_file, args.bigbed_file,
            args.sizes_file, args.output_file)


