import pandas as pd
import numpy as np
import itertools
import argparse

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-m', '--triangular-matrix', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-c', '--chrom', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-w', '--windows-file', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-o', '--output-format', default='npz', choices=['npz', 'threshold_npz', 'interaction_tsv', 'my5c_matrix'],
                    help='Format to convert to (options are "npz", "threshold_npz", or "interaction_tsv"')
parser.add_argument('-t', '--threshold-file', help='A file containing thresholds for calling confident interactions at each distance.')

def get_output_file(input_file, output_format):

    output_file = input_file.split('.')
    if output_format == 'npz':
        output_file[-1] = 'full_matrix.npz'
    elif output_format == 'threshold_npz':
        output_file[-1] = 'thresholded.npz'
    elif output_format == 'interaction_tsv':
        output_file[-1] = 'interactions.tsv'
    elif output_format == 'my5c_matrix':
        output_file[-1] = 'my5c.txt'
    else:
        raise Exception('Invalid choice of output format: {0}'.format(output_format))

    return '.'.join(output_file)

def open_triangular_matrix(filepath):
    with open(filepath) as in_data:
        arr = [[float(i) for i in line.split()] for line in in_data]
    N = len(arr[-1])
    full_array = np.zeros((N,N))
    lower_i = np.tril_indices_from(full_array)
    upper_i = np.triu_indices_from(full_array)
    full_array[:] = np.NAN
    full_array[lower_i] = list(itertools.chain(*arr))
    full_array[upper_i] = full_array.T[upper_i]
    full_array[full_array > 1.] = np.NAN
    
    return full_array

def open_windows(filepath, chrom):
    data = pd.read_csv(filepath, delim_whitespace=True, header=None)
    return np.array(data[data[0] == chrom])

def open_thresholds(filepath):
    return pd.read_csv(filepath, delim_whitespace=True, header=6).set_index('distance')

def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[:k], cols[-k:]
    elif k > 0:
        return rows[k:], cols[:-k]
    else:
        return rows, cols

def threshold_matrix_diagonals(matrix, thresholds):
    
    out_matr = np.zeros_like(matrix)
    
    for d in range(1, matrix.shape[0]+1):
        
        matrix_diag = matrix.diagonal(d).copy()
        
        try:
        
            thresh = np.array(thresholds.loc[d])
            
        except KeyError:
            
            thresh = np.array(thresholds.iloc[-1])
        
        below_thresh = matrix_diag < thresh
        
        matrix_diag[below_thresh] = 0.
        
        out_matr[kth_diag_indices(out_matr, d)] = matrix_diag
        out_matr[kth_diag_indices(out_matr, -d)] = matrix_diag
    
    return out_matr

def thresholded_to_df(matrix, chrom):

    _c = pd.DataFrame(matrix).unstack()
    
    _c = _c[_c > 0]
    _c =_c.reset_index()
    _c.columns = ['Pos_A', 'Pos_B', 'Pi']
    _c = _c[_c.Pos_B > _c.Pos_A]
    _c['dist'] = _c.Pos_B - _c.Pos_A
    _c['chrom'] = chrom
    return _c[['chrom', 'Pos_A', 'Pos_B', 'dist', 'Pi']]


if __name__ == '__main__':

    args = parser.parse_args()
    scores = open_triangular_matrix(args.triangular_matrix)
    windows = open_windows(args.windows_file, args.chrom)

    assert len(scores) == len(windows)

    output_file = get_output_file(args.triangular_matrix, args.output_format)

    if args.output_format == 'npz':
        np.savez_compressed(output_file, scores=scores, windows=windows)

    elif args.output_format == 'threshold_npz':
        if args.threshold_file is None:
            raise Exception('Threshold file not specified (use -t flag)')
        thresholds = open_thresholds(args.threshold_file)
        threshold_matrix = threshold_matrix_diagonals(scores, thresholds)
        np.savez_compressed(output_file, scores=threshold_matrix, windows=windows)

    elif args.output_format == 'interaction_tsv':
        if args.threshold_file is None:
            raise Exception('Threshold file not specified (use -t flag)')
        thresholds = open_thresholds(args.threshold_file)
        threshold_matrix = threshold_matrix_diagonals(scores, thresholds)
        interactions_df = thresholded_to_df(threshold_matrix, args.chrom)
        interactions_df.to_csv(output_file, sep='\t', index=False)


    elif args.output_format == 'my5c_matrix':
        names = [ '{}:{}-{}'.format(*i) for i in windows ]
        pd.DataFrame(scores, index=names, columns=names).to_csv(output_file, sep='\t', na_rep="NaN")

    else:
        raise Exception('Invalid choice of output format: {0}'.format(args.output_format))

