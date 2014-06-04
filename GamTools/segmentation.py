import numpy as np
import pandas as pd
import itertools
import os
from .cosegregation import Dprime


def open_segmentation(path_or_buffer):

    return pd.read_csv(path_or_buffer, 
                       index_col=[0,1,2], 
                       delim_whitespace=True)


def cosegregation_frequency(samples):
    """Take a table of n columns and return the co-segregation frequencies"""

    counts_shape = (2,) * samples.shape[0] 

    counts = np.zeros(counts_shape)

    for s in samples.T:

        counts[tuple(s)] += 1

    return counts


def get_index_combinations(regions):

    indexes = []
    start = 0

    for region in regions:
        indexes.append(range(start, start + len(region)))
        start = max(indexes[-1]) + 1

    return itertools.product(*indexes)


def get_cosegregation_freqs(*regions):

    if len(regions) == 1:

        regions = [ regions[0], regions[0] ]

    combinations = get_index_combinations(regions)

    full_data = np.concatenate(regions, axis=0)
    
    def get_frequency(indices):

        return cosegregation_frequency(full_data[indices, :])

    result = map(get_frequency, combinations)

    result_shape = tuple([ len(region) for region in regions ]) + (2, ) * len(regions)

    freqs = np.array(result).reshape(result_shape)

    return freqs


def map_freqs(func, fqs):
    old_shape = fqs.shape
    half_point = len(old_shape) / 2
    flat_shape = tuple([ np.prod(old_shape[:half_point]) ]) + old_shape[half_point:]
    return np.array(map(func, fqs.reshape(flat_shape))).reshape(old_shape[:half_point])


def index_from_interval(segmentation_data, interval):

    chrom, start, stop = interval

    window_in_region = np.logical_and(
                            np.logical_and(
                                segmentation_data.index.get_level_values('start') >= start,
                                segmentation_data.index.get_level_values('stop') <= stop),
                                segmentation_data.index.get_level_values('chrom') == chrom)

    covered_windows = np.nonzero(window_in_region)[0]
    start_index = covered_windows[0]
    stop_index = covered_windows[-1] + 1

    return start_index, stop_index


def parse_location_string(loc_string):

    chrom_fields = loc_string.split(':')

    chrom = chrom_fields[0]

    if len(chrom_fields) == 1:

        start, stop = 0, np.Inf

    else:

        pos_fields = chrom_fields[1].split('-')

        start, stop = map(int, pos_fields)

    return chrom, start, stop


def index_from_location_string(segmentation_data, loc_string):

    interval = parse_location_string(loc_string)

    return index_from_interval(segmentation_data, interval)


def region_from_location_string(segmentation_data, loc_string):

    ix_start, ix_stop = index_from_location_string(segmentation_data, loc_string)

    return segmentation_data.iloc[ix_start:ix_stop,]


def get_matrix(segmentation_data, *location_strings, **kwargs):

    def get_region(loc_string):

        return region_from_location_string(segmentation_data, loc_string)

    regions = map(get_region, location_strings)

    return get_matrix_from_regions(segmentation_data, *regions, **kwargs)


def get_matrix_from_regions(*regions, **kwargs):

    defaults = {'method' : Dprime,
               }

    defaults.update(kwargs)
    
    method = defaults['method']

    freqs = get_cosegregation_freqs(*regions)

    matrix = map_freqs(method, freqs)

    return matrix


def get_marginals(region):
    return region.sum(axis=1).astype(float) / region.count(axis=1).astype(float)


def map_sample_name_to_column(segmentation_data):

    name_mapping = { }

    for c in segmentation_data.columns:

        name_mapping[os.path.basename(c).split('.')[0]] = c  

    return name_mapping
