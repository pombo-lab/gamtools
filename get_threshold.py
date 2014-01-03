import numpy as np
import sys
try:
    import argparse
except ImportError:
    import backport_argparse as argparse
import itertools

parser = argparse.ArgumentParser(description='plot coverage over different window sizes for a list of bam files.')
parser.add_argument('input_file', metavar='INPUT_FILE', help='Input coverage file')
parser.add_argument('-t', '--threshold-file', metavar='THRESHOLD_FILE', help='File containing the threshold for each column')
parser.add_argument('--header-lines', metavar='NUMBER_OF_LINES', default=0, type=int, help='Skip this number of lines in the header')
parser.add_argument('--only-largest', action='store_true', help='Only output segmentation for the column with the largest number of reads')
args = parser.parse_args()

def get_data(fobj):
    data = []
    positions = []
    for line in itertools.islice(fobj, args.header_lines, None):
        fields = line.split()
        data.append(map(int,fields[3:]))
        positions.append([fields[0],'{0}-{1}'.format(*fields[1:3])])
    if args.only_largest:
        data = np.transpose(np.array(data))
        total_counts = map(max,data)
        biggest_count_index = total_counts.index(max(total_counts))
        return data[biggest_count_index], positions
    else:
        return data, positions

def get_threshold(fobj):

    reads_list = []
    threshold_list = []

    for line in fobj:
        reads, coverage, threshold = map(float,line.split())
        reads_list.append(reads)
        threshold_list.append(threshold)

    if args.only_largest:
        max_reads_index = reads_list.index(max(reads_list))

        return threshold_list[max_reads_index]
    else:
        return threshold_list

data, positions = get_data(open(args.input_file))
threshold_list = get_threshold(open(args.threshold_file))

def above_threshold(reads, threshold):
    if reads > threshold:
        return '1'
    else:
        return '0'

if args.only_largest:
    for i, reads in enumerate(data):
        above_threshold = '0'
        if reads > threshold:
            above_threshold = '1'
        fields = positions[i]
        fields.extend(above_threshold)
        print '\t'.join(fields)
else:
    sys.stdout.write('chrom\t')
    sys.stdout.write(open(args.input_file).readline())
    for i, read_list in enumerate(data):
        fields = positions[i]
        fields.extend(map(lambda p: above_threshold(read_list[p], threshold_list[p]), range(len(read_list))))
        print '\t'.join(fields)

