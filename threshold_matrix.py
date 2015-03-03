import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Apply a variable threshold to an npz matrix')
parser.add_argument('-n', '--npz-file', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-t', '--threshold-file', required=True, help='A file containing threshold values')
parser.add_argument('-s', '--skip', default=6, type=int, help='rows/columns at the end or beginning to skip')

args = parser.parse_args()

output_file = args.npz_file.split('.')
output_file = output_file[:-1] + ['thresholded'] + output_file[-1:]
output_file = '.'.join(output_file)

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

def load_npz(npz_path):

    windows = np.load(npz_path)['windows']
    scores = np.load(npz_path)['scores']

    return windows, scores

def load_threshold(threshold_path, skip=6):

    thresholds = pd.read_csv(threshold_path,
                delim_whitespace=True, header=skip).set_index('distance')

    return thresholds


windows, scores = load_npz(args.npz_file)
thresholds = load_threshold(args.threshold_file, args.skip)

thresholded_scores = threshold_matrix_diagonals(scores, thresholds)

np.savez_compressed(output_file, scores=thresholded_scores, windows=windows)

