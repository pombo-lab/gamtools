from .cosegregation_internal import cosegregation_2d, cosegregation_3d, linkage_2d, linkage_3d, dprime_2d
from . import segmentation, matrix
from .utils import format_genomic_distance
import numpy as np
import itertools
import time
import warnings
import sys


class InvalidDataError(Exception):
    """Exception to raise if segmentation data contains anything other than 0s and 1s"""
    pass


def regions_are_valid(regions):

    invalid_regions = [(np.array(region) > 1).any() for region in regions]

    if any(invalid_regions):
        return False
    else:
        return True


def prepare_regions(regions):

    if not regions_are_valid(regions):
        raise InvalidDataError('Region contains integers greater than 1')

    if len(regions) == 1:
        regions = regions * 2

    regions = [np.array(r.astype(int)) for r in regions]

    return regions


def _cosegregation_frequency(samples):
    """Take a table of n columns and return the co-segregation frequencies"""

    counts_shape = (2,) * samples.shape[0]

    counts = np.zeros(counts_shape)

    for s in samples.T:

        counts[tuple(s)] += 1

    return counts


def get_index_combinations(regions):

    indexes = []
    start = 0

    assert len(regions) > 1

    for region in regions:

        indexes.append(range(start, start + len(region)))
        start = max(indexes[-1]) + 1

    return itertools.product(*indexes)


def cosegregation_nd(*regions):

    combinations = get_index_combinations(regions)

    full_data = np.concatenate(regions, axis=0).astype(int)

    def get_frequency(indices):

        return _cosegregation_frequency(full_data[indices, :])

    result = map(get_frequency, combinations)

    result_shape = tuple([ len(region) for region in regions ]) + (2, ) * len(regions)

    freqs = np.array(result).reshape(result_shape)

    return freqs


def get_cosegregation_from_regions(*regions):

    regions = prepare_regions(regions)

    if len(regions) == 2:
        coseg_func = cosegregation_2d
    elif len(regions) == 3:
        coseg_func = cosegregation_3d
    else:
        coseg_func = cosegregation_nd

    return coseg_func(*regions)


def get_cosesgregation(segmentation_data, *location_strings):

    regions = [segmentation.region_from_location_string(segmentation_data, l) for l in location_strings]

    return get_cosegregation_from_regions(*regions)


def map_freqs(func, fqs):
    old_shape = fqs.shape
    half_point = len(old_shape) / 2
    flat_shape = tuple([ np.prod(old_shape[:half_point]) ]) + old_shape[half_point:]
    return np.array(map(func, fqs.reshape(flat_shape))).reshape(old_shape[:half_point])


def get_linkage_from_regions(*regions):

    regions = prepare_regions(regions)

    if len(regions) == 2:
        linkage_func = linkage_2d
    elif len(regions) == 3:
        warnings.warn('3D linkage is calculated following Hastings - Genetics (1984) 106:153-164. Please check and make sure this is what you want.')
        linkage_func = linkage_3d
    else:
        raise NotImplementedError('No linkage function is defined for {} regions'.format(len(regions)))

    return linkage_func(*regions)


def get_linkage(segmentation_data, *location_strings):

    regions = [segmentation.region_from_location_string(segmentation_data, l) for l in location_strings]

    return get_linkage_from_regions(*regions)


def get_dprime_from_regions(*regions):

    regions = prepare_regions(regions)

    if len(regions) == 2:
        dprime_func = dprime_2d
    else:
        raise NotImplementedError('There is currently no implementation of normalized linkage disequilibrium for more than 2 dimensions')

    return dprime_func(*regions)


def get_dprime(segmentation_data, *location_strings):

    regions = [segmentation.region_from_location_string(segmentation_data, l) for l in location_strings]

    return get_dprime_from_regions(*regions)

matrix_types = {
    'dprime': get_dprime_from_regions,
    'linkage': get_linkage_from_regions,
    'cosegregation': get_cosegregation_from_regions,
}

def get_regions_and_windows(segmentation_data, location_strings):

    if len(location_strings) == 1:
        location_strings = location_strings * 2

    regions = [segmentation.region_from_location_string(segmentation_data, location_str)
               for location_str in location_strings]

    windows = [np.array(list(region.index)) for region in regions]

    return regions, windows

def matrix_and_windows_from_segmentation_file(
    segmentation_file, location_strings, matrix_type='dprime'):

    segmentation_data = segmentation.open_segmentation(segmentation_file)

    regions, windows = get_regions_and_windows(segmentation_data, location_strings)
    matrix_func = matrix_types[matrix_type]
    contact_matrix = matrix_func(*regions)

    return contact_matrix, windows


def create_and_save_contact_matrix(segmentation_file, location_strings,
                                   output_file, output_format,
                                   matrix_type='dprime'):

    print 'starting calculation for {}'.format(' x '.join(location_strings))
    start_time = time.clock()

    contact_matrix, windows = matrix_and_windows_from_segmentation_file(
        segmentation_file, location_strings, matrix_type)

    size_string = ' x '.join([str(s) for s in contact_matrix.shape])
    print 'region size is: {}'.format(size_string),
    print 'Calculation took {0}s'.format(time.clock() - start_time)
    print 'Saving matrix to file {}'.format(output_file)

    output_func = matrix.output_formats[output_format]

    output_func(windows, contact_matrix, output_file)

    print 'Done!'

def get_output_file(segmentation_file, location_strings, matrix_type, output_format):

    segmentation_base = '.'.join(segmentation_file.split('.')[:-1])
    locations = [segmentation.parse_location_string(loc_string)
                 for loc_string in location_strings]
    formatted_locations = [(chrom, format_genomic_distance(start), format_genomic_distance(stop))
                           for chrom, start, stop in locations]
    formatted_locations = ['{}_{}-{}'.format(*loc) for loc in formatted_locations]
    region_string = '_x_'.join(formatted_locations)

    return '{}.{}_{}.{}'.format(segmentation_base, region_string, matrix_type, output_format)

def matrix_from_args(args):

    if args.output_format is None:
        if len(args.regions) > 2:
            args.output_format = 'npz'
        elif args.output_file is None:
            args.output_format = 'txt.gz'
        elif args.output_file == '-':
            args.output_format = 'txt'
        else:
            args.output_format = matrix.detect_file_type(args.output_file)

    if args.output_file is None:
        args.output_file = get_output_file(
            args.segmentation_file, args.regions,
            args.matrix_type, args.output_format)

    elif args.output_file == '-':
        args.output_file = sys.stdout

    create_and_save_contact_matrix(args.segmentation_file, args.regions,
                                   args.output_file, args.output_format,
                                   args.matrix_type)
