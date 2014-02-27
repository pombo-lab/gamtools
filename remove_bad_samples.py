import pandas as pd
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Remove columns in segmentation file that fail QC.')
parser.add_argument('-s','--segmentation-file', metavar='SEGMENTATION_FILE', required=True, help='A file containing the segmentation of all samples')
parser.add_argument('-n','--sample-names', metavar='SAMPLE_NAME', required=True, nargs='*', help='Names of the samples to remove')

args = parser.parse_args()

data = pd.read_csv(args.segmentation_file, sep='\s*')

index = { os.path.basename(c).split('.')[0] : c for c in data.columns }

bad_samples = [ index[b] for b in args.sample_names ]

subset = data.drop(bad_samples, 1)

subset.to_csv(sys.stdout, index=False, sep='\t')

