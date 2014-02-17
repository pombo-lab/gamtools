import argparse
import pandas as pd
from itertools import groupby

parser = argparse.ArgumentParser(description='Takes in a segmentation multibam file and returns proportion of gaps equal to min(gap) for each column')
parser.add_argument('segmentation', metavar='SEGMENTATION_MULTIBAM', help='Path to a segmentation multibam')

def proportion_smallest_gap(block_list):
    gaps = []
    for k, g in groupby(block_list):
        if k:
            gaps.extend([0] * ( len(list(g)) - 1))
        else:
            gaps.append(len(list(g)))
            
    return float(len([gap for gap in gaps if gap == min(gaps)])) / len(gaps)

if __name__ == '__main__':
    args = parser.parse_args()

    data = pd.read_csv(args.segmentation, delim_whitespace=True)
    for c in data.columns[2:]:
        print c, proportion_smallest_gap(data[c])
