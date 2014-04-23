from GamTools import GamExperiment
import itertools
import numpy as np
import argparse
from scipy.stats import pearsonr

parser = argparse.ArgumentParser(description='Calculate PCC between different chromosomes stored as npz files.')
parser.add_argument('-s1','--segmentation1', metavar='FIRST_SEGMENTATION', required=True, help='First segmentation multibam')
parser.add_argument('-s2','--segmentation2', metavar='SECOND_SEGMENTATION', required=True, help='Second segmentation multibam')

if __name__ == '__main__':

    args = parser.parse_args()

    exp1 = GamExperiment.from_multibam_no_cache(args.segmentation1)
    exp2 = GamExperiment.from_multibam_no_cache(args.segmentation2)

    chroms = [ c for c in set(map(lambda x:x[0], exp1.experimental_data.windows)) if c[-7:] not in ['_random', 'chrM', 'chrY'] ]

    # Make chromosome iterators
    gen1 = ( exp1.distances(c).flatten() for c in chroms )
    gen2 = ( exp2.distances(c).flatten() for c in chroms )

    # Chain chromosome iterators
    chain1 = np.array(list(itertools.chain.from_iterable(gen1)))
    chain2 = np.array(list(itertools.chain.from_iterable(gen2)))

    na1 = np.ma.getmask(np.ma.masked_invalid(chain1))
    na2 = np.ma.getmask(np.ma.masked_invalid(chain2))
    neitherna = np.logical_not(np.logical_or(na1, na2))

    print 'Correlation is: {0}'.format(pearsonr(np.compress(neitherna, chain1), np.compress(neitherna, chain2))[0])
