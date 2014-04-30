import pandas as pd
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Remove columns in segmentation file that fail QC.')
parser.add_argument('-s','--segmentation-file', metavar='SEGMENTATION_FILE', required=True, help='A file containing the segmentation of all samples')
parser.add_argument('-i','--interactions-file', metavar='INTERACTION_FILE', required=True, help='A file containing the interactions to calculate')

args = parser.parse_args()

data = pd.read_csv(args.segmentation_file,
                   delim_whitespace=True, index_col=[0,1,2]).transpose()

def get_chrom_start_stop_bins(chrom):

    chrom_bins = map(lambda t: t[0], data.columns)

    start_bin = chrom_bins.index(chrom)
    stop_from_end = chrom_bins[::-1].index(chrom)
    stop_bin = len(chrom_bins) - 1 - stop_from_end

    return start_bin, stop_bin


def count_frequency(samples):
    """Take a table of n columns and return the co-segregation frequencies"""

    counts_shape = (2,) * samples.shape[1] 
    
    counts = np.zeros(counts_shape)

    for s in samples:

        counts[tuple(s)] += 1

    return counts


def get_transpositions(array):
    axes = range(len(array.shape))
    for i in axes:
        yield tuple(axes[i:] + axes[:i])


def get_independent_probs(n):
    
    total = n.sum()
    f = n / float(total)
    
    ind = []
    for t in get_transpositions(f):
        probs = [ f.transpose(t)[1,...].sum(), f.transpose(t)[0,...].sum() ]
        ind.append(probs)
    return np.array(ind), f


def D(n):
    probs,f = get_independent_probs(n)
    expected = probs.prod(axis=0)[0]
    observed = f.flatten()[-1]
    return observed - expected

def corr(n):
    d =  D(n)
    probs,f = get_independent_probs(n)
    return d / np.power(probs.prod(), 1./len(n.shape))


def line_indices(line):
    fields = line.split()
    chrom = fields[0]
    cstart, cstop = get_chrom_start_stop_bins(chrom)
    indices = map(lambda x: cstart+int(x),fields[1:])
    
    assert all([i < cstop for i in indices])
    
    return indices
    
def indices_corr(indices):
    
    columns = np.array(data.iloc[:, indices ])
    freqs = count_frequency(columns)
    
    return corr(freqs), freqs

def line_corr(line):
    
    indices = line_indices(line)
    correlation, frequencies = indices_corr(indices)
    return (correlation,) + tuple(get_independent_probs(frequencies)[0][:,0])
    
with open(args.interactions_file, 'r') as interactions:
    for line in interactions:
        correlation = line_corr(line)
        original_fields = line.strip().split()
        correlation_fields = list(correlation)
        correlation_fields.append(max(correlation_fields[1:]))
        print '\t'.join(map(str,original_fields + correlation_fields))
