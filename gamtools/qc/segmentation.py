import pandas as pd
from .. import segmentation
import os
import itertools

chroms = ['chr{0}'.format(c) for c in range(1,20)]

def proportion_with_neighbours(block_list):
    no_with_neighbours = 0
    for k, g in itertools.groupby(block_list):
        g = list(g)
        if k and len(g) > 1:
            no_with_neighbours += len(g)

    try:
        return float(no_with_neighbours) / sum(block_list)
    except ZeroDivisionError:
        return 0.

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

def get_segmentation_stats(input_segmentation, output_file):

    segmentation_df = segmentation.open_segmentation(input_segmentation)

    stats_df = get_df_stats(segmentation_df)

    stats_df.to_csv(output_file, sep='\t', index=False)

def get_segmentation_stats_doit(dependencies, targets):

    assert len(dependencies) == 1
    assert len(targets) == 1

    get_segmentation_stats(list(dependencies)[0], list(targets)[0])
