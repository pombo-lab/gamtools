import numpy as np
import pandas as pd
import itertools
import os
import glob
from .cosegregation import Dprime
from .coseg_3 import cosegregation_frequency_3, coseg_explicit_2_cy, coseg_explicit_3_cy
#from EIYBrowse.filetypes.my5c_folder import My5CFolder, My5cFile

class InvalidChromError(Exception):
    """Exception to be raised when an invalid chromosome is specified"""
    pass

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


def fast_cosegregation_frequency(samples):
    """Take a table of 2 columns and return the co-segregation frequencies"""

    coseg = np.product(samples, axis=0).sum()
    marga, margb = np.sum(samples, axis=1)

    return np.array([[len(samples[0]) - marga - margb + coseg,
                      margb - coseg],[
                      marga - coseg, 
                      coseg]])

def get_index_combinations(regions):

    indexes = []
    start = 0

    assert len(regions) > 1

    for region in regions:

        indexes.append(range(start, start + len(region)))
        start = max(indexes[-1]) + 1

    return itertools.product(*indexes)

def get_cosegregation_freqs(*regions):

    if len(regions) == 1:
        regions = regions * 2

    if len(regions) == 2:
        coseg_func = fast_cosegregation_frequency
    elif len(regions) == 3:
        coseg_func = cosegregation_frequency_3
    else:
        coseg_func = cosegregation_frequency

    combinations = get_index_combinations(regions)

    full_data = np.concatenate(regions, axis=0).astype(int)
    
    def get_frequency(indices):

        return coseg_func(full_data[indices, :])

    result = map(get_frequency, combinations)

    result_shape = tuple([ len(region) for region in regions ]) + (2, ) * len(regions)

    freqs = np.array(result).reshape(result_shape)

    return freqs

def coseg_explicit_2(seg1, seg2):

    counts = np.zeros((2,2))

    for i in range(len(seg1)):

        counts[seg1[i], seg2[i]] += 1

    return counts

def get_cosegregation_freqs_explicit_2(region_1, region_2):

    result = []
    res_append = result.append

    for xi, yi in itertools.product(range(len(region_1)), range(len(region_2))):
        res_append(coseg_explicit_2_cy(region_1[xi], region_2[yi]))

    result_shape = (len(region_1), len(region_2), 2, 2)

    freqs = np.array(result).reshape(result_shape)

    return freqs


def get_cosegregation_freqs_explicit_3(region_1, region_2, region_3):

    result = []
    res_append = result.append

    for xi, yi, zi in itertools.product(xrange(len(region_1)),
                                        xrange(len(region_2)),
                                        xrange(len(region_3))):

        res_append(coseg_explicit_3_cy(region_1[xi],
                                       region_2[yi],
                                       region_3[zi]))

    result_shape = (len(region_1), len(region_2), len(region_3),
                    2, 2, 2)

    freqs = np.array(result).reshape(result_shape)

    return freqs

def map_freqs(func, fqs):
    old_shape = fqs.shape
    half_point = len(old_shape) / 2
    flat_shape = tuple([ np.prod(old_shape[:half_point]) ]) + old_shape[half_point:]
    return np.array(map(func, fqs.reshape(flat_shape))).reshape(old_shape[:half_point])


def index_from_interval(segmentation_data, interval):

    chrom, start, stop = interval

    if not start < stop:
        raise ValueError('Interval start {0} larger than interval end {1}'.format(*interval))

    window_in_region = np.logical_and(
                            np.logical_and(
                                segmentation_data.index.get_level_values('stop') > start,
                                segmentation_data.index.get_level_values('start') < stop),
                                segmentation_data.index.get_level_values('chrom') == chrom)

    covered_windows = np.nonzero(window_in_region)[0]

    if not len(covered_windows):
        if not chrom in segmentation_data.index.levels[0]:
            raise InvalidChromError('{0} not found in the list of windows'.format(chrom))

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
    def __init__(self, segmentation_path, method=Dprime):
        super(GamSegmentationFile, self).__init__()

        self.data = open_segmentation(segmentation_path)
        self.method = method
    
    def interactions(self, feature):
        
        interval = feature.chrom, feature.start, feature.stop
        ix_start, ix_stop = index_from_interval(self.data, interval)
        region = self.data.iloc[ix_start:ix_stop,]
                
        return get_matrix_from_regions(region, method=self.method), feature

class TooManyFilesError(Exception):
    pass

class NoFilesError(Exception):
    pass

"""
class GamNpzFolder(My5CFolder):
    def __init__(self, folder_path):
        
        self.folder_path = folder_path
        self.file_class = GamNpzFile
        
    def find_chrom_file(self, chrom):
        
        search_string = '*{0}[_\.]*npz'.format(chrom)
        
        found_files = glob.glob(os.path.join(self.folder_path, search_string))
        
        if len(found_files) > 1:
            raise TooManyFilesError('folder containing npzs must have only one npz per chromosome. Found {0}'.format(found_files))
        elif not found_files:
            raise NoFilesError('No npz file found matching {0} with search string "{1}"'.format(chrom, search_string))
            
        return found_files[0]
    
class GamNpzFile(My5cFile):

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
    """
