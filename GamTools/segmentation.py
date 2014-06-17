import numpy as np
import pandas as pd
import itertools
import os
import glob
from .cosegregation import Dprime
import pybedtools


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

        regions = regions * 2

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
                                segmentation_data.index.get_level_values('stop') >= start,
                                segmentation_data.index.get_level_values('start') <= stop),
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

    return get_matrix_from_regions(*regions, **kwargs)


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

class GamSegmentationFile(object):
    """Panel for displaying a continuous signal (e.g. ChIP-seq) across a genomic region"""
    def __init__(self, segmentation_path):
        super(GamSegmentationFile, self).__init__()

        self.data = open_segmentation(segmentation_path)
    
    def interactions(self, feature):
        
        interval = feature.chrom, feature.start, feature.stop
        ix_start, ix_stop = index_from_interval(self.data, interval)
        region = self.data.iloc[ix_start:ix_stop,]
                
        return get_matrix_from_regions(region), feature

class TooManyFilesError(Exception):
    pass

class GamNpzFolder(object):
    def __init__(self, folder_path):
        
        self.folder_path = folder_path
        
    def find_chrom_file(self, chrom):
        
        search_string = '*.{0}.npz'.format(chrom)
        
        found_files = glob.glob(os.path.join(self.folder_path, search_string))
        
        if len(found_files) > 1:
            raise TooManyFilesError('folder containing npzs must have only one npz per chromosome')
            
        return found_files[0]

    def get_npz_file(self, chrom):

        npz_path = self.find_chrom_file(chrom)

        return GamNpzFile(npz_path)
    
    def interactions(self, feature):
        
        npz_file = self.get_npz_file(feature.chrom)
        
        return npz_file.get_interactions(feature)
    
class GamNpzFile(object):
    def __init__(self, file_path):
        
        data = np.load(file_path)
        self.interactions = data['scores']
        self.windows = self.format_windows(data['windows'])
        
    def format_windows(self, windows):
        
        def format_window(window):
            
            chrom = window[0]
            start, stop = map(int, window[1:])
            return chrom, start, stop
        
        return pd.MultiIndex.from_tuples(map(format_window, windows),
                                         names=['chrom', 'start', 'stop'])
    
    def index_from_interval(self, feature):

        window_in_region = np.logical_and(
                                np.logical_and(
                                    self.windows.get_level_values('start') >= feature.start,
                                    self.windows.get_level_values('stop') <= feature.stop),
                                    self.windows.get_level_values('chrom') == feature.chrom)

        covered_windows = np.nonzero(window_in_region)[0]
        start_index = covered_windows[0]
        stop_index = covered_windows[-1] + 1 

        return start_index, stop_index

    def get_interactions(self, feature):
        
        start, stop = self.index_from_interval(feature)
        
        return self.interactions[start:stop, start:stop], self.indices_to_interval(start, stop)
    
    def indices_to_interval(self, start, stop):
        
        start_window = self.windows[start]
        stop_window = self.windows[stop - 1]
        return pybedtools.Interval(start_window[0], start_window[1], stop_window[2])
