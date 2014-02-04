import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Takes in a segmentation multibam file and returns coverage stats for each column')
parser.add_argument('segmentation', metavar='SEGMENTATION_MULTIBAM', help='Path to a segmentation multibam')

def get_stats(multibam_handle):
    data = pd.read_csv(multibam_handle, sep='\s*')
    summed = data.groupby('chrom').sum()
    print 'Sample', 'Total_positive_windows', 'Chromosomes_above_mean'
    for col in summed.columns:
        print col, data[col].sum(), len(summed[summed[col] > summed[col].mean()])

if __name__ == '__main__':
    args = parser.parse_args()

    with open(args.segmentation, 'r') as multibam_handle:
        get_stats(multibam_handle)
