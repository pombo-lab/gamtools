import sys
import argparse
from GamTools import segmentation

parser = argparse.ArgumentParser(description='Remove columns in segmentation file that fail QC.')
parser.add_argument('-s','--segmentation-file', metavar='SEGMENTATION_FILE', required=True, help='A file containing the segmentation of all samples')
parser.add_argument('-g','--keep-good', action='store_true', help='Keep only the listed samples (default: discard listed samples)')
parser.add_argument('-n','--sample-names', metavar='SAMPLE_NAME', required=True, nargs='*', help='Names of the samples to remove')

args = parser.parse_args()

if __name__ == '__main__':

    data = segmentation.open_segmentation(args.segmentation_file)

    index = segmentation.map_sample_name_to_column(data)

    bad_samples = [ index[b] for b in args.sample_names ]

    if args.keep_good:
        subset = data[bad_samples]
    else:
        subset = data.drop(bad_samples, 1)


    subset.to_csv(sys.stdout, index=True, sep='\t')

