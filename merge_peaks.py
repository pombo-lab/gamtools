import argparse
import pandas as pd
import sys
import os

parser = argparse.ArgumentParser(description='plot coverage over different window sizes for a list of bam files.')
parser.add_argument('-w', '--windows', metavar='WINDOWS_FILE', required=True, help='Input coverage file')
parser.add_argument('input_files', metavar='INPUT_FILE', nargs='+', help='Input coverage file')

args = parser.parse_args()

def sub_name(windows, new_name):
    w = list(windows.columns)
    w[w.index('pos')] = new_name
    windows.columns = w

windows = pd.read_csv(args.windows, delim_whitespace=True,
                header=None, names=['chrom', 'start', 'stop'], index_col=[0,1,2])

for bed in args.input_files:
    bed_data = pd.read_csv(bed, delim_whitespace=True,
                header=None, names=['chrom', 'start', 'stop', 'pos'], index_col=[0,1,2])
    windows = pd.merge(windows, bed_data, left_index=True, right_index=True, how='left')
    sub_name(windows, os.path.basename(bed).split('.')[0])

windows = windows.fillna(0)
windows.astype(int).to_csv(sys.stdout, index=True, sep='\t')
