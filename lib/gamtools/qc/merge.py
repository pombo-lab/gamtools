"""
===================
The qc.merge module
===================

The qc.merge module contains functions for merging separate qc stats files
into a single table.

"""
import pandas as pd

def merge_stats(input_stats_files, output_merged_file):
    """
    Merge a list of dataframes together based on their index columns.

    :param list input_stats_files: String paths to input dataframes.
    :param str output_merged_file: Path to save output dataframe.
    """

    first_file = input_stats_files[0]

    base = pd.read_csv(first_file, delim_whitespace=True, index_col=0)

    for stats_file_path in input_stats_files[1:]:

        stats_file = pd.read_csv(stats_file_path, delim_whitespace=True, index_col=0)
        base = pd.merge(base, stats_file, left_index=True, right_index=True)

    base.to_csv(output_merged_file, sep='\t')

def merge_stats_from_doit(dependencies, targets):
    """
    Wrapper function to call merge_stats from argparse.
    """

    assert len(targets) == 1
    merge_stats(list(dependencies), targets[0])
