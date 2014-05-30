import numpy as np
import pandas as pd
import itertools


def open_segmentation(path_or_buffer):

    return pd.read_csv(path_or_buffer, 
                       index_col=[0,1,2], 
                       delim_whitespace=True)


def cosegregation_frequency(samples):
    """Take a table of n columns and return the co-segregation frequencies"""

    counts_shape = (2,) * samples.shape[1] 
    
    counts = np.zeros(counts_shape)

    for s in samples:

        counts[tuple(s)] += 1

    return counts


def get_index_combinations(regions):

    indexes = []
    start = 0

    for region in regions:
        indexes.append(range(start, start + len(region.columns)))
        start = max(indexes[-1]) + 1

    return itertools.product(*indexes)


def get_cosegregation_freqs(*regions):

    if len(regions) == 1:

        regions = [ regions[0], regions[0] ]

    combinations = get_index_combinations(regions)

    full_data = np.concatenate(regions, axis=1)
    
    def get_frequency(indices):
                
        return cosegregation_frequency(full_data[:, indices])

    result = map(get_frequency, combinations)

    result_shape = tuple([ len(region.columns) for region in regions ]) + (2, ) * len(regions)

    freqs = np.array(result).reshape(result_shape)

    return freqs


def map_freqs(func, fqs):
    old_shape = fqs.shape
    half_point = len(old_shape) / 2
    flat_shape = tuple([ np.prod(old_shape[:half_point]) ]) + old_shape[half_point:]
    return np.array(map(func, fqs.reshape(flat_shape))).reshape(old_shape[:half_point])

def get_marginals(region):
    return region.sum(axis=1).astype(float) / region.count(axis=1).astype(float)
