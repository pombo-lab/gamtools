import numpy as np
import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Circularly permute every column of a segmentation file')
parser.add_argument('segmentation_file', help='A segmentation file to use as input')

def permute_data(input_data):
    no_windows, no_samples = input_data.shape
    permutation = input_data.copy()
    for i in range(no_samples):
        offset = np.random.randint(no_windows)
        permutation.iloc[:,i] = list(pd.concat([input_data.iloc[offset:,i], input_data.iloc[:offset,i]]))
    return permutation

if __name__ == '__main__':

    args = parser.parse_args()
    
    segmentation = pd.read_csv(args.segmentation_file,
                          delim_whitespace=True, index_col=[0,1,2])

    permutation = permute_data(segmentation)

    permutation.to_csv(sys.stdout, sep='\t', header=True, index=True)
