import numpy as np
import pandas as pd
import sys, os
from .segregation import parse_location_string
import itertools
import argparse

def get_name_strings(windows):

    return ['{}:{}-{}'.format(*i) for i in windows]

def windows_from_name_strings(name_strings):

    return [parse_location_string(name) for name in name_strings]

def read_npz(filepath):

    handle = np.load(filepath)
    contact_matrix = handle['scores']

    try:
        windows = [handle['windows_{}'.format(i)] for i in range(contact_matrix.ndim)]
    except KeyError:
        windows = [handle['windows'], handle['windows']]

    return windows, contact_matrix

def read_txt(filepath, compression='infer'):

    contact_matrix = pd.read_csv(filepath, compression=compression, sep='\t', index_col=0)

    windows_0 = windows_from_name_strings(contact_matrix.index)
    windows_1 = windows_from_name_strings(contact_matrix.columns)

    return (windows_0, windows_1), contact_matrix

def read_zipped_txt(filepath):

    return read_txt(filepath, compression='gzip')

def read_windows(filepath, chrom):

    data = pd.read_csv(filepath, delim_whitespace=True, header=None)

    windows = np.array(data[data[0] == chrom])

    return [windows, windows]

def read_triangular(filepath):

    with open(filepath) as in_data:
        arr = [[float(i) for i in line.split()] for line in in_data]
    N = len(arr[-1])
    contact_matrix = np.zeros((N,N))
    lower_i = np.tril_indices_from(contact_matrix)
    upper_i = np.triu_indices_from(contact_matrix)
    contact_matrix[:] = np.NAN
    contact_matrix[lower_i] = list(itertools.chain(*arr))
    contact_matrix[upper_i] = contact_matrix.T[upper_i]
    contact_matrix[contact_matrix > 1.] = np.NAN

    return None, contact_matrix


input_formats = {
        'npz': read_npz,
        'txt': read_txt,
        'txt.gz': read_txt,
        'triangular': read_triangular,
        # TODO: Convert from interactions csv back into a matrix
        #'csv': '',
}

def write_txt(windows, contact_matrix, output_file):

    if contact_matrix.ndim != 2:
        raise NotImplementedError('Plain text output is only supported for 2 dimensional matrices. Please try saving as an npz file.')

    names_0, names_1 = [get_name_strings(axis_windows) for axis_windows in windows]

    pd.DataFrame(contact_matrix, index=names_0, columns=names_1).to_csv(output_file, sep='\t', na_rep="NaN")

def write_zipped_txt(windows, contact_matrix, output_file):

    import gzip
    with gzip.open(output_file, 'wb', compresslevel=5) as zipped_output:
        write_txt(windows, contact_matrix, zipped_output)

def write_npz(windows, contact_matrix, output_file):

    window_dict = {'windows_{}'.format(i):win for i,win in enumerate(windows)}
    np.savez_compressed(output_file, scores=contact_matrix, **window_dict)

def write_csv(windows, contact_matrix, output_file):

        interactions_df = pd.DataFrame(contact_matrix).unstack()
        interactions_df = interactions_df[interactions_df > 0]
        interactions_df =interactions_df.reset_index()
        interactions_df.columns = ['Pos_A', 'Pos_B', 'interaction']
        interactions_df = interactions_df[interactions_df.Pos_B > interactions_df.Pos_A]
        interactions_df['dist'] = interactions_df.Pos_B - interactions_df.Pos_A
        interactions_df['chrom'] = windows[0][0][0]
        output_cols = ['chrom', 'Pos_A', 'Pos_B', 'dist', 'interaction']
        interactions_df[output_cols].to_csv(output_file, sep='\t', index=False)

def write_zipped_csv(windows, contact_matrix, output_file):

    import gzip
    with gzip.open(output_file, 'wb', compresslevel=5) as zipped_output:
        write_csv(windows, contact_matrix, zipped_output)

output_formats = {
        'npz': write_npz,
        'txt': write_txt,
        'txt.gz': write_zipped_txt,
        'csv': write_csv,
        'csv.gz': write_zipped_csv,
}

def detect_file_type(file_name):

    if file_name == '-':
        return 'txt'

    file_name = os.path.basename(file_name)
    file_parts = file_name.split('.')

    if len(file_parts) == 1:
        raise ValueError('Could not determine file format, file {} has no extension'.format(file_name))

    file_ext = file_parts[-1]

    if file_ext == 'gz':
        file_ext = '.'.join(file_parts[-2:])

    if file_ext in output_formats.keys():
        return file_ext

    raise TypeError('Extension "{}" not recognized'.format(file_ext))


def check_windows(contact_matrix, windows):

    windows_sizes = [len(win) for win in windows]
    for i in range(contact_matrix.ndim):
        if contact_matrix.shape[i] != windows_sizes[i]:
            raise ValueError(
                'Contact matrix size ({}) does not match the number '
                'of windows supplied ({}). Please check you have '
                'specified the correct region.'.format(
                    ' x '.join([str(s) for s in contact_matrix.shape]),
                    ' x '.join([str(s) for s in windows_sizes])))


def read_thresholds(thresholds_file):
    return pd.read_csv(thresholds_file,
                       delim_whitespace=True, header=6).set_index('distance')


def kth_diag_indices(a, k):
        rows, cols = np.diag_indices_from(a)
        if k < 0:
            return rows[:k], cols[-k:]
        elif k > 0:
            return rows[k:], cols[:-k]
        else:
            return rows, cols


def apply_threshold(contact_matrix, thresholds):

    out_matr = np.zeros_like(contact_matrix)

    for d in range(1, contact_matrix.shape[0]+1):
        contact_matrix_diag = contact_matrix.diagonal(d).copy()

        try:
            thresh = np.array(thresholds.loc[d])
        except KeyError:
            thresh = np.array(thresholds.iloc[-1])

        below_thresh = contact_matrix_diag < thresh
        contact_matrix_diag[below_thresh] = 0.
        out_matr[kth_diag_indices(out_matr, d)] = contact_matrix_diag
        out_matr[kth_diag_indices(out_matr, -d)] = contact_matrix_diag

    return out_matr


def convert(input_file, input_format,
            output_file, output_format,
            windows=None, thresholds=None):

    _windows, contact_matrix = input_formats[input_format](input_file)

    if input_format == 'triangular':
        if windows is None:
            raise argparse.ArgumentError(
                None,
                'A windows file must be specified (-w/--windows-file)'
                ' when converting triangular matrices')

        _windows = windows

    check_windows(contact_matrix, _windows)

    if thresholds is not None:
        contact_matrix = apply_threshold(contact_matrix, thresholds)

    output_formats[output_format](_windows, contact_matrix, output_file)

def convert_from_args(args):

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
