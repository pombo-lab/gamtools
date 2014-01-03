import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser(description='Extract from a segmentation file only the positive windows for a particular sample')
parser.add_argument('-s','--segmentation-file', metavar='SEGMENTATION_FILE', required=True, help='A file containing the segmentation of all samples')
parser.add_argument('-n','--sample-name', metavar='SAMPLE_NAME', required=True, help='Name of the sample to extract')

args = parser.parse_args()

data = pd.read_csv(args.segmentation_file, sep='\s*')

subset = data.loc[data[args.sample_name] > 0,['chrom','window',args.sample_name]]

subset.to_csv(sys.stdout, header=False, index=False, sep='\t')
