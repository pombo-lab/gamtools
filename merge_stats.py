import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Takes in a segmentation multibam file and returns proportion of gaps equal to min(gap) for each column')
parser.add_argument('stats_files', metavar='STATS_FILE', nargs='+', help='Path to a segmentation multibam')

if __name__ == '__main__':
    args = parser.parse_args()

    first_file = args.stats_files.pop()

    base = pd.read_csv(first_file, delim_whitespace=True, index_col=0)

    for stats_file_path in args.stats_files:

        stats_file = pd.read_csv(stats_file_path, delim_whitespace=True, index_col=0)
        base = pd.merge(base, stats_file, left_index=True, right_index=True)

    base.to_csv(sys.stdout, sep='\t')
