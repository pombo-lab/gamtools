import pandas as pd
import itertools
import numpy as np
import argparse
import logging
import os
import sys

parser = argparse.ArgumentParser(description='Calculate z-score for a chromosome stored as My5C files.')
parser.add_argument('-m','--my5c-path', metavar='MY5C_PATH', required=True, help='An input My5C file.')
parser.add_argument('--debug',
    help='Print lots of debugging statements',
    action="store_const",dest="loglevel",const=logging.DEBUG,
    default=logging.WARNING
)
parser.add_argument('--verbose',
    help='Be verbose',
    action="store_const",dest="loglevel",const=logging.INFO
)

args = parser.parse_args()

chrom_data = np.array(pd.read_csv(args.my5c_path, sep='\t', index_col=0))
chrom_data[chrom_data == 0.] = np.NAN
chrom_data = np.log10(chrom_data)

zscore_file = args.my5c_path.split('.')
zscore_file.insert(-1,"zscore")
zscore_file = '.'.join(zscore_file)

if os.path.exists(zscore_file):
    sys.exit('File exists already, not overwriting')

def get_mean_std(chrom_data):
    indices = itertools.product(xrange(chrom_data.shape[0]), xrange(chrom_data.shape[1]))
    table_data = [ (p, q, chrom_data[p,q]) for p,q in indices ] 
    df = pd.DataFrame.from_records(table_data, columns=['x', 'y', 'correlation'])
    df['dist'] = np.abs(df.x - df.y)
    dists = df.groupby('dist')['correlation']
    means = dists.mean()
    stdevs = dists.std()
    return means, stdevs

means, stdevs = get_mean_std(chrom_data)

def get_z_score(x, dist):
    return (x - means[dist]) / stdevs[dist]

zs = np.array( [ get_z_score(chrom_data[p,q], np.abs(p - q)) for p,q in 
                        itertools.product(xrange(chrom_data.shape[0]),
                                          xrange(chrom_data.shape[1])) ] )

np.savez_compressed(zscore_file, scores=zs.reshape(chrom_data.shape))
