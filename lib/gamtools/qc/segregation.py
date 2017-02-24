"""
=========================
The qc.segregation module
=========================

The qc.segregation module contains functions for calculating quality control
statistics from segregation tables.

"""

import os
import itertools

import pandas as pd

from .. import segregation

def proportion_with_neighbours(block_list):
    """
    Calculate the percentage of positive windows that have positive neighbours

    :param list block_list: list of positive windows.
    :returns: Percentage positive windows with positive neighbours.
    """

    no_with_neighbours = 0
    for key, group in itertools.groupby(block_list):
        group = list(group)
        if key and len(group) > 1:
            no_with_neighbours += len(group)

    try:
        return float(no_with_neighbours) / sum(block_list)
    except ZeroDivisionError:
        return 0.

def extract_sample_name(path):
    """
    Get the sample name from the path
    """

    return os.path.basename(path).split('.')[0]

def get_df_stats(segregation_df):
    """
    Generate a table of summary statistics for each NP in a segregation table.

    Statistics calculated are as follows:

    - *Proportion_with_neighbours*: Percentage of positive windows with a
      positive neighbour
    - *Genome_coverage*: Percentage of all windows which are positive
    - *Positive_chromosomes*: Number of chromosomes with greater than
      mean percentage of positive windows.

    :param segregation_df: Input :ref:`~segregation table`
    :returns: :class:`~pandas.DataFrame` of statistics for each NP.
    """

    prop_neighb = segregation_df.apply(proportion_with_neighbours, axis=0)
    prop_neighb.name = 'Proportion_with_neighbours'

    genome_cover = segregation_df.mean()
    genome_cover.name = 'Genome_coverage'

    positive_chroms = (segregation_df.groupby(level=0).mean() >
                       segregation_df.groupby(level=0).mean().mean()).sum()
    positive_chroms.name = 'Positive_chromosomes'

    stats_df = pd.concat([genome_cover, positive_chroms, prop_neighb], axis=1).reset_index()

    stats_df['Sample'] = stats_df['index'].apply(extract_sample_name)

    return stats_df[['Sample', 'Genome_coverage', 'Positive_chromosomes',
                     'Proportion_with_neighbours']]

def get_segregation_stats(input_segregation, output_file):
    """
    Given the path to a segregation table file, open the file and save the
    summary statistics to output_file.

    :param str input_segregation: Path to input segregation file.
    :param str output_file: Path to save the output file.
    """

    segregation_df = segregation.open_segregation(input_segregation)

    stats_df = get_df_stats(segregation_df)

    stats_df.to_csv(output_file, sep='\t', index=False)

def get_segregation_stats_doit(dependencies, targets):
    """Wrapper function to call get_segregation_stats from argparse"""

    assert len(dependencies) == 1
    assert len(targets) == 1

    get_segregation_stats(list(dependencies)[0], list(targets)[0])
