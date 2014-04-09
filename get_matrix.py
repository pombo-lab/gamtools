from mpi4py import MPI
import itertools
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Map GAM-seq reads and create BigWigs and fastQC files')
parser.add_argument('-o','--output_file', metavar='OUTPUT_FILE', required=True, help='Write frequencies to this file')
parser.add_argument('-c','--chromosome', metavar='CHROMOSOME', required=True, help='Specific chromosome to calculate matrices for')
parser.add_argument('segmentation_file', help='A segmentation file to use as input')

args = parser.parse_args()

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

data = pd.read_csv(args.segmentation_file,
                   delim_whitespace=True, index_col=[0,1,2]).transpose()

def get_chrom_start_stop_bins(chrom):

    chrom_bins = map(lambda t: t[0], data.columns)

    start_bin = chrom_bins.index(chrom)
    stop_from_end = chrom_bins[::-1].index(chrom)
    stop_bin = len(chrom_bins) - 1 - stop_from_end

    return start_bin, stop_bin

loc1_start, loc1_stop = get_chrom_start_stop_bins(args.chromosome)
loc2_start, loc2_stop = loc1_start, loc1_stop

data_1 = np.array(data.iloc[:, loc1_start:loc1_stop])
len_1 = len(data_1[0])
data_2 = np.array(data.iloc[:, loc2_start:loc2_stop])
len_2 = len(data_2[0])
full_data = np.concatenate((data_1, data_2), axis=1)

def count_frequency(samples):
    """Take a table of two columnds and return [[ no_both_present, no_1_only],[ no_2_only, no_neither_present]]"""

    counts = np.array([[0, 0], [0, 0]])

    for s in samples:

        counts[s[0]][s[1]] += 1

    return counts

def get_frequency(i):

    p, q = i
    return count_frequency(full_data[:, [p, q]])

rows = range(len_1)
 
cols = range(len_1, len_1 + len_2)

my_rows = itertools.ifilter(lambda x: x % size == rank, rows)

product = []
for r in my_rows:
    row = itertools.product([r], cols)
    result = map(get_frequency, list(row))
    product.append(result)

full = comm.gather(product, root=0)

if rank == 0:

    zipped = itertools.izip_longest(*full)
    no_None = [ list(itertools.ifilter(lambda p: p is not None, pair)) for pair in zipped ]
    freqs = np.array(list(itertools.chain(*no_None)))
    np.savez_compressed(args.output_file, freqs=freqs)
