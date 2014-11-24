import argparse
import pandas as pd
from GamTools import segmentation
import os
import sys
import itertools

parser = argparse.ArgumentParser(description='Takes in a segmentation multibam file and returns coverage stats for each column')
parser.add_argument('segmentation', metavar='SEGMENTATION_MULTIBAM', help='Path to a segmentation multibam')

chroms = ['chr{0}'.format(c) for c in range(1,20)]

def proportion_with_neighbours(block_list):
    no_with_neighbours = 0
    for k, g in itertools.groupby(block_list):
        g = list(g)
        if k and len(g) > 1:
            no_with_neighbours += len(g)
            
    return float(no_with_neighbours) / sum(block_list)

def extract_sample_name(path):
    return os.path.basename(path).split('.')[0]

def get_df_stats(segmentation_df):

    prop_neighb = segmentation_df.apply(proportion_with_neighbours, axis=0)
    prop_neighb.name = 'Proportion_with_neighbours'

    genome_cover = segmentation_df.mean()
    genome_cover.name = 'Genome_coverage'

    positive_chroms = (segmentation_df.groupby(level=0).mean() > segmentation_df.groupby(level=0).mean().mean()).sum()
    positive_chroms.name = 'Positive_chromosomes'

    stats_df = pd.concat([genome_cover, positive_chroms, prop_neighb], axis=1).reset_index()

    stats_df['Sample'] = stats_df['index'].apply(extract_sample_name)

    return stats_df[['Sample', 'Genome_coverage', 'Positive_chromosomes',
            'Proportion_with_neighbours']]

def print_file_stats(multibam_handle):

    segmentation_df = segmentation.open_segmentation(multibam_handle)

    stats_df = get_df_stats(segmentation_df)

    stats_df.to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
    args = parser.parse_args()

    with open(args.segmentation, 'r') as multibam_handle:
        print_file_stats(multibam_handle)
