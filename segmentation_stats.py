import argparse
from GamTools.segmentation import open_segmentation
from itertools import groupby
import numpy as np
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Takes in a segmentation multibam file and returns proportion of gaps equal to min(gap) for each column')
parser.add_argument('segmentation', metavar='SEGMENTATION_MULTIBAM', help='Path to a segmentation multibam')

def proportion_smallest_gap(block_list):
    no_with_neighbours = 0
    for k, g in groupby(block_list):
        g = list(g)
        if k and len(g) > 1:
            no_with_neighbours += len(g)
            
    return float(no_with_neighbours) / len(block_list)

if __name__ == '__main__':
    args = parser.parse_args()

    data = open_segmentation(args.segmentation)
    result_df = pd.DataFrame(np.array(data.sum(axis=0).astype(float) / data.count(axis=0)), 
                             columns=['total_pos_windows'], index=data.columns)

    result_df.index.name = 'Sample'

    result_df['pct_adjacent'] = [proportion_smallest_gap(np.array(data[c])) for c in data.columns ]

    autosomal = ['chr{0}'.format(c) for c in range(0,20)]
    data_by_chrom = data[[c in autosomal for c in data.index.get_level_values(0)]].groupby(level=0)
    chrom_coverages = data_by_chrom.sum().astype(float) / data_by_chrom.count()
    result_df['chrom_above_mean'] = (chrom_coverages.T > list(chrom_coverages.mean())).sum(axis=1)

    result_df.to_csv(sys.stdout, sep='\t')
