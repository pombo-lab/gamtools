"""
=================
The matrix module
=================

The matrix module contains functions for reading and writing proximity
matrices in various formats.

.. _proximity_matrices:

Proximity matrices
==================

One of the most common ways of displaying and processing information about
chromatin folding is as a 2D proximity matrix, where the x-axis represents the
position along a contiguous genomic region and the y-axis represents the
position along another contiguous genomic region (which can be the same region
as on the x-axis). Each point within the matrix gives the inferred distance
between location x and location y, with higher values generally indicating a
shorter distance. For example, this type of matrix is the main output of a Hi-C
experiment (Lieberman-Aiden et al. 2009, Science).


.. _matrix_formats:

Proximity matrix formats
========================

npz
---

Npz format is the native compressed numpy format. It takes up less disk space
and is faster to read and write than the other formats, but it is not human
readable and is not compatible with other programming languages like R.

If you want to work with data that involves interactions between more
than two partners (i.e. three for four-way interactions), that data can be
represented as an n-dimensional matrix (tensor). Npz format is the only
format supported by GAMtools that can read or write n-dimensional matrices.

txt
---

Text file format is the simplest matrix format. It is a tab-delimited text file
where the columns and rows are the columns and rows of the matrix, plus one header
row and one column of row names. Text files are human readable but take up a lot
of disk space and are slow to read and write. GAMtools can also read and write
to gzipped text files by using the txt.gz format (these take up much less disk
space).

triangular
----------

Triangular matrix files are similar to text files but they only give one half
of a symmetrical matrix, and they do not contain column/row names. They are
the format used by SLICE.

csv
---

csv files are text files used as input to the enrichment module. The format
can only be used for saving intrachromosomal data. The first column gives
the chromosome name, the second and third columns give the column and row
of the matrix (as 0 based indices, not by using their genomic locations)
the fourth column gives the distance between the x and y co-ordinate and the
fifth column gives the interaction strength. See :func:`write_csv` for an
example.

png
---

Specifying the png output format will save the proximity matrix to a 2D
heatmap image.


"""

import sys
import os
import itertools
import argparse

import numpy as np
import pandas as pd

from . import segregation
from .utils import DelayedImportError

try:
    from matplotlib import pyplot as plt
except ImportError:
    plt = DelayedImportError(
        'Saving a matrix as an image requires matplotlib to be installed. '
        'Try to install it by running "pip install matplotlib"')


def get_name_strings(windows):
    """Format a list of tuples representing genomic windows (chrom, start, stop)
    as a list of UCSC-like location strings "chrom:start-stop"

    :param list windows: List of tuples in the format (chromosome, start_coord, stop_coord)
    :returns: List of UCSC formatted location strings
    """

    return ['{}:{}-{}'.format(*i) for i in windows]


def windows_from_name_strings(name_strings):
    """Convert a list of UCSC-like location strings to a list of window tuples

    :param list name_strings: List of UCSC-like location strings
    :returns: List of window tuples
    """

    return [segregation.parse_location_string(name) for name in name_strings]


def read_npz(filepath):
    """Open an npz file containing a proximity matrix

    :param str filepath: Path to the npz file
    :returns: List of lists giving genomic locations for each bin on each axis, and a\
            :class:`numpy array <numpy.ndarray>` proximity matrix.
    """

    handle = np.load(filepath)
    proximity_matrix = handle['scores']

    try:
        windows = [
            handle[
                'windows_{}'.format(i)] for i in range(
                    proximity_matrix.ndim)]
    except KeyError:
        windows = [handle['windows'], handle['windows']]

    return windows, proximity_matrix


def read_txt(filepath, compression='infer'):
    """Open an txt file containing a proximity matrix

    :param str filepath: Path to the txt file
    :param str compression: One of {'infer', 'gzip', 'bz2', None}, \
    default 'infer'. For on-the-fly decompression of on-disk data. If \
    'infer', then use gzip or bz2 if filepath is a string \
    ending in '.gz' or '.bz2', respectively, and no decompression \
    otherwise. Set to None for no decompression.
    :returns: (list of genomic locations for x-axis, list of \
    genomic locations for y-axis), \
            :class:`numpy array <numpy.ndarray>` proximity matrix.
    """

    proximity_matrix = pd.read_csv(
        filepath,
        compression=compression,
        sep='\t',
        index_col=0)

    windows_0 = windows_from_name_strings(proximity_matrix.index)
    windows_1 = windows_from_name_strings(proximity_matrix.columns)

    return (windows_0, windows_1), proximity_matrix.values


def read_windows(filepath, chrom):
    """Retrieve the genomic locations (windows) covering a given chromosome
    from a text file.

    This function is used when opening Pi matrices output from SLICE, as
    these matrices do not have embedded window locations.

    :param str filepath: Path to the text file
    :param str chrom: Chromosome to fetch the windows for.
    :returns: (list of genomic locations for x-axis, list of \
    genomic locations for y-axis)
    """

    data = pd.read_csv(filepath, delim_whitespace=True, header=None)

    windows = np.array(data[data[0] == chrom])

    return [windows, windows]


def read_triangular(filepath):
    """Open Pi matrix output from SLICE.

    All matrix opening functions return first the
    genomic windows corresponding to the axes of the
    proximity matrix, then the proximity matrix itself.
    Since SLICE output matrices do not embed the genomic
    locations of the windows, the first return value is
    None.

    :param str filepath: Path to the SLICE output file
    :returns: (None, SLICE Pi matrix)
    """

    with open(filepath) as in_data:
        arr = [[float(i) for i in line.split()] for line in in_data]
    size = len(arr[-1])
    proximity_matrix = np.zeros((size, size))
    lower_i = np.tril_indices_from(proximity_matrix)
    upper_i = np.triu_indices_from(proximity_matrix)
    proximity_matrix[:] = np.NAN
    proximity_matrix[lower_i] = list(itertools.chain(*arr))
    proximity_matrix[upper_i] = proximity_matrix.T[upper_i] # pylint: disable=E1136
    proximity_matrix[proximity_matrix > 1.] = np.NAN

    return None, proximity_matrix


INPUT_FORMATS = {
    'npz': read_npz,
    'txt': read_txt,
    'txt.gz': read_txt,
    'triangular': read_triangular,
    # TODO: Convert from interactions csv back into a matrix
    #'csv': '',
}


def write_txt(windows, proximity_matrix, output_file):
    """Write a proximity matrix to a txt file.

    :param tuple windows: (list of x-axis windows, list of y-axis windows)
    :param proximity_matrix: Input proximity matrix.
    :type proximity_matrix: :class:`numpy array <numpy.ndarray>`
    :param str filepath: Path to save matrix file.

    >>> my_matrix = np.array([[10, 0, 5],
    ...                       [0, 10, 3],
    ...                       [5, 3, 10]])
    >>> windows = [('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 30)]
    >>> matrix.write_txt([windows, windows], my_matrix, sys.stdout)
        chr1:0-10       chr1:10-20      chr1:20-30
        chr1:0-10       10      0       5
        chr1:10-20      0       10      3
        chr1:20-30      5       3       10

    """

    if proximity_matrix.ndim != 2:
        raise NotImplementedError(
            'Plain text output is only supported for 2 dimensional matrices. '
            'Please try saving as an npz file.')

    names_0, names_1 = [get_name_strings(
        axis_windows) for axis_windows in windows]

    pd.DataFrame(
        proximity_matrix, index=names_0, columns=names_1).to_csv(
            output_file, sep='\t', na_rep="NaN")


def write_npz(windows, proximity_matrix, output_file):
    """Write a proximity matrix to an npz file.

    npz files are a compressed numpy-specific format, meaning
    they take up less disk space, but cannot be easily opened
    by other programming languages (e.g. R). For more information
    see :func:`numpy.savez_compressed`.

    :param tuple windows: (list of x-axis windows, list of y-axis windows)
    :param proximity_matrix: Input proximity matrix.
    :type proximity_matrix: :class:`numpy array <numpy.ndarray>`
    :param str filepath: Path to save matrix file.
    """

    window_dict = {
        'windows_{}'.format(i): win for i,
        win in enumerate(windows)}
    np.savez_compressed(output_file, scores=proximity_matrix, **window_dict)


def write_csv(windows, proximity_matrix, output_file):
    """Write a proximity matrix to a csv file.

    csv file outputs a table giving the index of the windows on the x- and
    y-axis, the interaction score (if the score is greater than 0), the
    chromosome, and the distance between the two windows. It is only
    appropriate for intra-chromosomal proximity matrices because the
    matrix is assumed to be symmetrical, and duplicated information
    is discarded.

    :param tuple windows: (list of x-axis windows, list of y-axis windows)
    :param proximity_matrix: Input proximity matrix.
    :type proximity_matrix: :class:`numpy array <numpy.ndarray>`
    :param str filepath: Path to save csv file.

    >>> my_matrix = np.array([[10, 0, 5],
    ...                       [0, 10, 3],
    ...                       [5, 3, 10]])
    >>> windows = [('chr1', 0, 10), ('chr1', 10, 20), ('chr1', 20, 30)]
    >>> matrix.write_csv([windows, windows], my_matrix, sys.stdout)
    chrom   Pos_A   Pos_B   dist    interaction
    chr1    0       2       2       5
    chr1    1       2       1       3
    """

    interactions_df = pd.DataFrame(proximity_matrix).unstack()
    interactions_df = interactions_df[interactions_df > 0]
    interactions_df = interactions_df.reset_index()
    interactions_df.columns = ['Pos_A', 'Pos_B', 'interaction']
    interactions_df = interactions_df[
        interactions_df.Pos_B > interactions_df.Pos_A]
    interactions_df['dist'] = interactions_df.Pos_B - interactions_df.Pos_A
    interactions_df['chrom'] = windows[0][0][0]
    output_cols = ['chrom', 'Pos_A', 'Pos_B', 'dist', 'interaction']
    interactions_df[output_cols].to_csv(output_file, sep='\t', index=False)


def write_png(windows, proximity_matrix, output_file): #pylint: disable=unused-argument
    """Write a proximity matrix to a .png image file.

    :param tuple windows: (list of x-axis windows, list of y-axis windows)
    :param proximity_matrix: Input proximity matrix.
    :type proximity_matrix: :class:`numpy array <numpy.ndarray>`
    :param str filepath: Path to save image file.
    """

    plt.figure(figsize=(7, 7))
    plt.imshow(proximity_matrix, interpolation='none', cmap='jet')
    plt.axis('off')
    plt.savefig(output_file)

OUTPUT_FORMATS = {
    'npz': write_npz,
    'txt': write_txt,
    'txt.gz': write_txt,
    'csv': write_csv,
    'csv.gz': write_csv,
    'png': write_png,
}

SUPPORTED_FORMATS = list(set(list(INPUT_FORMATS.keys()) + list(OUTPUT_FORMATS.keys())))


def detect_file_type(file_name):
    """Given the path to a matrix file, determine the file type
    based on the extension.

    :param str file_name: Path to the matrix file
    :returns: File type string
    :raises TypeError: If file type does not match one of the keys \
            in either the output_formats or input_formats dictionaries.

    >>> matrix.detect_file_type('my_file.txt')
    'txt'
    >>> matrix.detect_file_type('my_file.png')
    'png'
    >>> matrix.detect_file_type('my_file.data')
    TypeError: Extension "doc" not recognized
    """

    if file_name == '-':
        return 'txt'

    file_name = os.path.basename(file_name)
    file_parts = file_name.split('.')

    if len(file_parts) == 1:
        raise ValueError(
            'Could not determine file format, file {} has no extension'.format(file_name))

    file_ext = file_parts[-1]

    if file_ext == 'gz':
        file_ext = '.'.join(file_parts[-2:])

    if file_ext in SUPPORTED_FORMATS:
        return file_ext

    raise TypeError('Extension "{}" not recognized'.format(file_ext))


def read_file(file_name, file_type=None):
    """Open a matrix file, guessing the format based on file extension.

    :param str file_name: Path to matrix file.
    :returns: (list of genomic locations for x-axis, list of \
    genomic locations for y-axis), \
            :class:`numpy array <numpy.ndarray>` proximity matrix.
    """

    if file_type is None:
        file_type = detect_file_type(file_name)
    read_func = INPUT_FORMATS[file_type]
    return read_func(file_name)


def check_windows(proximity_matrix, windows):
    """Check that a list of axis windows matches the dimensions of a proximity matrix.

    :param proximity_matrix: Input proximity matrix.
    :type proximity_matrix: :class:`numpy array <numpy.ndarray>`
    :param tuple windows: (list of x-axis windows, list of y-axis windows)
    :raises ValueError: If the number of windows doesn't match the proximity matrix dimensions.
    """

    windows_sizes = [len(win) for win in windows]
    for i in range(proximity_matrix.ndim):
        if proximity_matrix.shape[i] != windows_sizes[i]:
            raise ValueError(
                'Contact matrix size ({}) does not match the number '
                'of windows supplied ({}). Please check you have '
                'specified the correct region.'.format(
                    ' x '.join([str(s) for s in proximity_matrix.shape]),
                    ' x '.join([str(s) for s in windows_sizes])))


def region_from_locations(matrix_tuple, location_string1, location_string2=None):
    """Retrieve a sub-region of an interaction matrix based on a location string.

    :param matrix_tuple: Tuple of the form (x_windows_list, y_windows_list), matrix
    :param str location_string1: UCSC format :ref:`location string <location_string>`
    :param str location_string2: UCSC format :ref:`location string <location_string>`, \
            if None, location_string2 is set equal to location_string1 (i.e. the function \
            returns a symmetrical sub-matrix)
    :returns: Tuple of the form (x_windows_list, y_windows_list), matrix
    """

    (windows1, windows2), chr_mat = matrix_tuple

    windows1_df = pd.DataFrame(
        windows1, columns=['chrom', 'start', 'stop']).set_index(['chrom', 'start', 'stop'])
    windows2_df = pd.DataFrame(
        windows2, columns=['chrom', 'start', 'stop']).set_index(['chrom', 'start', 'stop'])

    if location_string2 is None:
        location_string2 = location_string1

    i1_start, i1_stop = segregation.index_from_location_string(windows1_df, location_string1)
    i2_start, i2_stop = segregation.index_from_location_string(windows2_df, location_string2)

    subregion_w1 = windows1_df.iloc[i1_start:i1_stop]
    subregion_w2 = windows2_df.iloc[i2_start:i2_stop]

    submatrix = chr_mat[i1_start:i1_stop, i2_start:i2_stop]

    return (subregion_w1, subregion_w2), submatrix


def open_region_from_locations(matrix_path,
                               location_string1, location_string2=None,
                               file_type=None):
    """Open an interaction matrix and retrieve a sub-region based on a location string.

    :param matrix_path: Path to an interaction matrix
    :param str location_string1: UCSC format :ref:`location string <location_string>`
    :param str location_string2: UCSC format :ref:`location string <location_string>`, \
            if None, location_string2 is set equal to location_string1 (i.e. the function \
            returns a symmetrical sub-matrix)
    :returns: Tuple of the form (x_windows_list, y_windows_list), matrix
    """

    matrix_tuple = read_file(matrix_path, file_type)

    return region_from_locations(matrix_tuple, location_string1, location_string2)


def read_thresholds(thresholds_file):
    """Read a file containing interaction thresholds"""

    return pd.read_csv(thresholds_file,
                       delim_whitespace=True, comment='#').set_index('distance')


def kth_diag_indices(array, diag_k):
    """Return a tuple of indices for retrieving the k'th diagonal
    of matrix a.

    :param array: Input matrix.
    :type array: :class:`numpy array <numpy.ndarray>`
    :param int diag_k: Diagonal to index. 0 is the centre, 1 is the first \
            diagonal below the centre, -1 is the first diagonal above \
            the centre.

    >>> my_matrix = np.array([[ 0, -1, -2, -3],
    ...                       [ 1,  0, -1, -2],
    ...                       [ 2,  1,  0, -1],
    ...                       [ 3,  2,  1,  0]])
    >>> matrix.kth_diag_indices(my_matrix, 1)
    (array([1, 2, 3]), array([0, 1, 2]))
    >>> my_matrix[
    ...    matrix.kth_diag_indices(my_matrix, 1)
    ...  ]
    array([1, 1, 1])
    >>> my_matrix[
    ...    matrix.kth_diag_indices(my_matrix, -2)
    ...  ]
    array([-2, -2])
    """

    rows, cols = np.diag_indices_from(array)
    if diag_k < 0:
        return rows[:diag_k], cols[-diag_k:]

    if diag_k > 0:
        return rows[diag_k:], cols[:-diag_k]

    return rows, cols


def apply_threshold(proximity_matrix, thresholds):
    """Discard values in a proximity matrix which are below a
    threshold.

    Thresholds are applied separately to each diagonal (i.e.
    each distance between windows). Thus, the first threshold
    in the list is applied to the first diagonal (i.e. to all
    adjacent windows). If the number of thresholds is less
    than the number of diagonals, the last threshold is repeated.

    :param proximity_matrix: Input proximity matrix.
    :type proximity_matrix: :class:`numpy array <numpy.ndarray>`
    :param thresholds: Values to use as minimum thresholds.
    :type thresholds: :class:`pandas.DataFrame`
    """

    out_matr = np.zeros_like(proximity_matrix)

    for diagonal in range(1, proximity_matrix.shape[0] + 1):
        proximity_matrix_diag = proximity_matrix.diagonal(diagonal).copy()

        try:
            thresh = np.array(thresholds.loc[diagonal])
        except KeyError:
            thresh = np.array(thresholds.iloc[-1])

        below_thresh = proximity_matrix_diag < thresh
        proximity_matrix_diag[below_thresh] = 0.
        out_matr[kth_diag_indices(out_matr, diagonal)] = proximity_matrix_diag
        out_matr[kth_diag_indices(out_matr, -diagonal)] = proximity_matrix_diag

    return out_matr


def convert(input_file, input_format, #pylint: disable=too-many-arguments
            output_file, output_format,
            windows=None, thresholds=None):
    """Convert a proximity matrix to a different file format.

    :param str input_file: Path to input proximity matrix
    :param str input_format: File format of input matrix
    :param str output_file: Path to output proximity matrix
    :param str output_format: File format to convert matrix to
    :param list windows: List of genomic windows (only required \
            when converting from SLICE output matrices)
    :param thresholds: Table of interaction thresholds (index \
            column is the distance to apply the threshold \
            given in number of bins)
    :type thresholds: :class:`pandas.DataFrame`
    """

    _windows, proximity_matrix = INPUT_FORMATS[input_format](input_file)

    if input_format == 'triangular':
        if windows is None:
            raise argparse.ArgumentError(
                None,
                'A windows file must be specified (-w/--windows-file)'
                ' when converting triangular matrices')

        _windows = windows

    check_windows(proximity_matrix, _windows)

    if thresholds is not None:
        proximity_matrix = apply_threshold(proximity_matrix, thresholds)

    OUTPUT_FORMATS[output_format](_windows, proximity_matrix, output_file)


def convert_from_args(args):
    """Wrapper function to call convert from argparse"""

    if args.input_format is None:
        args.input_format = detect_file_type(args.input_file)

    if args.output_format is None:
        args.output_format = detect_file_type(args.output_file)

    if args.input_file == '-':
        args.input_file = sys.stdin

    if args.output_file == '-':
        args.output_file = sys.stdout

    if args.windows_file is not None:
        if args.region is None:
            raise argparse.ArgumentError(
                None,
                'A region must be specified (-r/--region)'
                ' when converting triangular matrices')

        windows = read_windows(args.windows_file, args.region)
    else:
        windows = None

    if args.thresholds_file is not None:
        thresholds = read_thresholds(args.thresholds_file)
    else:
        thresholds = None

    convert(args.input_file, args.input_format,
            args.output_file, args.output_format,
            windows, thresholds)
