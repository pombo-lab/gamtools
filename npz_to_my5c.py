import numpy as np
import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Convert gam heatmap in npz format to txt format.')
parser.add_argument('-n', '--npz_frequencies_file', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')

args = parser.parse_args()

def open_npz(fp):
    handle = np.load(fp)
    return handle['windows'], handle['scores']

windows, data = open_npz(args.npz_frequencies_file)
names = [ '{}:{}-{}'.format(*i) for i in windows ]
pd.DataFrame(data, index=names, columns=names).to_csv(sys.stdout, sep='\t', na_rep="NaN")
