import argparse
import pandas as pd
from itertools import groupby

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

    data = pd.read_csv(args.segmentation, delim_whitespace=True)
    print data
    for c in data.columns[2:]:
        print c, proportion_smallest_gap(data[c])
