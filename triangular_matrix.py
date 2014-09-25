import pandas as pd
import numpy as np
import itertools
import argparse

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-t', '--triangular-matrix', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-c', '--chrom', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-w', '--windows-file', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')

args = parser.parse_args()

output_file = args.triangular_matrix.split('.')
output_file[-1] = 'full_matrix.npz'
output_file = '.'.join(output_file)

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

scores = open_triangular_matrix(args.triangular_matrix)
windows = open_windows(args.windows_file, args.chrom)

assert len(scores) == len(windows)

np.savez_compressed(output_file, scores=scores, windows=windows)
