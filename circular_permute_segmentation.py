import numpy as np
import argparse
import pandas as pd
from GamTools.segmentation import open_segmentation
import sys

parser = argparse.ArgumentParser(description='Circularly permute every column of a segmentation file')
parser.add_argument('segmentation_file', help='A segmentation file to use as input')

def permute_data(input_data):

    # Only permute positions that were mapped at least once
    mappable = input_data.sum(axis=1).astype(bool)
    no_windows, no_samples = input_data[mappable].shape

    # Make a copy of the original data
    permutation = input_data.copy()

    # Loop over columns
    for i in range(no_samples):

        # Choose a random position to break the segmentation in two,
        # Swap the two chunks around and write them to the copied df
        offset = np.random.randint(no_windows)
        
        new_start = input_data[mappable].iloc[offset:,i]
        new_end = input_data[mappable].iloc[:offset,i]
        new_col = list(pd.concat([new_start, new_end]))
        
        permutation.ix[mappable,i] = new_col
    return permutation

if __name__ == '__main__':

    args = parser.parse_args()
    
    segmentation = open_segmentation(args.segmentation_file)

    # Sex chromosomes, unjoined contigs and mtDNA behave weirdly,
    # so don't permute them into the autosomal regions.

    autosomes = ['chr{0}'.format(c) for c in range(1,20)]
    
    segmentation = segmentation.loc[autosomes]

    permutation = permute_data(segmentation)

    permutation.to_csv(sys.stdout, sep='\t', header=True, index=True)
