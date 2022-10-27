"""
========================
The cosegregation module
========================

Measuring proximity in GAM datasets
===================================

The fundamental principle behind GAM is that two regions which are in close
proximity in the nucleus should be frequently found in the same thin
nuclear slice - i.e. that more proximal genomic regions co-segregate more
frequently than less proximal genomic regions. The cosegregation module
provides functions that calculate the co-segregation of genomic regions,
thereby allowing the relative nuclear proximity of different regions
to be inferred.

Co-segregation, the number of times that location x and location y are found in
the same :ref:`nuclear profile <NPs>`, is the simplest measure of proximity but
there are alternative approaches. Other ways of generating a :ref:`proximity
matrix <proximity_matrices>` generally attempt to normalize for the
differential detection of the two locations. For example, if locations x and y
are detected in 100 :ref:`NPs` and location z is detected in only 20, the
co-segregation of x and y will likely be higher than that between x and z even
if their respective nuclear proximities are the same. We have previously
reported that the approach which best accounts for such detection effects is
the :func:`normalized linkage disequilibrium <get_dprime>` (or Dprime).
However, as of GAMtools v2.0 the recommended metric (and the new default for
all matrix calculations) is :func"`normalized pointwise mutual information \
        <get_npmi>`.


"""
from __future__ import print_function
import itertools
import time
import warnings
import sys

import numpy as np

from .cosegregation_internal import cosegregation_2d, cosegregation_3d, \
        linkage_2d, linkage_3d, dprime_2d
from .npmi import npmi_2d
from . import segregation, matrix
from .utils import format_genomic_distance


class InvalidDataError(Exception):
    """
    Exception raised if segregation data contains anything other than 0s and 1s.
    """


def regions_are_valid(regions):
    """
    Test whether any regions contain values other than the integers 0 and 1

    :param list regions: List of :ref:`regions <regions>` to check.
    :returns: Returns False if any regions contain invalid data, otherwise True.
    """

    allowed_values = set([0, 1])
    region_unique_values = [
        set(np.unique(np.array(region).ravel())) for region in regions]
    invalid_regions = [not region_values.issubset(
        allowed_values) for region_values in region_unique_values]

    if any(invalid_regions):
        return False

    return True


def prepare_regions(regions):
    """Checks and formats a list of regions

    Takes a list of :ref:`regions <regions>`, checks if they are valid
    (see :func:`regions_are_valid`) and converts them to integer arrays.
    If there is only one region in the list, cosegregation should be calculated
    for that region against itself, so returns a list containing the same
    region twice

    :param list regions: List of :ref:`regions <regions>` to check.
    :returns: list of integer :class:`numpy arrays <numpy.ndarray>`.
    """

    if not regions_are_valid(regions):
        raise InvalidDataError('Region contains integers greater than 1')

    if len(regions) == 1:
        regions = regions * 2

    regions = [np.array(r.astype(int)) for r in regions]

    return regions


def cosegregation_frequency_ndim(loci):
    """Take a table of n columns and return the n-dimensional co-segregation frequencies.

    This function accepts an integer :class:`numpy array <numpy.ndarray>` of size n by x, where
    n is the number of genomic windows (loci) and x is the number of samples. It calculates
    the number of times from 0 to n of the windows are all found in the same tube.

    This is a pure python function, so it is relatively slow, but it can calculate
    the co-segregation of any number of windows.

    :param loci: Array giving the segregation of n loci across x samples.
    :type loci: :class:`~numpy.ndarray`
    :returns: An n-dimensional contingency table giving the cosegregation frequencies of all \
    possible combinations of the n loci.
    """

    counts_shape = (2,) * loci.shape[0]

    counts = np.zeros(counts_shape)

    for sample in loci.T:

        counts[tuple(sample)] += 1

    return counts


def get_index_combinations(regions):
    """Returns all possible index combinations for a list of regions.

    Takes a list of :ref:`regions <regions>` and returns a generator
    that will iterate through indices for all combinations of columns
    in the three regions. Column indices for each region are consecutive,
    for example, if the first region has 20 columns, the last column of
    region 1 will have the index 19 (since indices are 0-based) and the
    first column of region 2 will have the index 20.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: generator that yields tuples of integer indices.
    """

    indexes = []
    start = 0

    assert len(regions) > 1

    for region in regions:

        indexes.append(list(range(start, start + len(region))))
        start = max(indexes[-1]) + 1

    return itertools.product(*indexes)


def cosegregation_nd(*regions):
    """Get the full co-segregation table for n regions.

    Takes a list of :ref:`regions <regions>` and calculates the
    full cosegregation frequency matrix for all possible combinations
    of loci from the n regions. For example, if three regions are
    given, position (x,y,z) of the output matrix is a three-dimensional
    contingency table, giving the number of times locus x cosegregates
    with loci y and z.

    This function has a pure python inner and outer loop, so it is
    very slow, but it has the advantage of being applicable in an
    unlimited number of dimensions. Optimized functions for
    obtaining the co-segregation of two or three regions can
    be found in the cosegregation_optimized module.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the co-segregation \
            of all possible combinations of windows within the different regions.
    """

    combinations = get_index_combinations(regions)

    full_data = np.concatenate(regions, axis=0).astype(int)

    def get_frequency(indices):
        """
        Internal function to get the cosegregation_frequency for
        a given combination of columns from full_data.
        """

        return cosegregation_frequency_ndim(full_data[indices, :])

    result = list(map(get_frequency, combinations))

    result_shape = tuple(len(region)
                         for region in regions) + (2, ) * len(regions)

    freqs = np.array(result).reshape(result_shape)

    return freqs


def get_cosegregation_from_regions(*regions):
    """Get the full co-segregation table for n regions, using optimized
    functions where available.

    This is a wrapper which determines the correct co-segregation
    function to call based on the number of regions. Where there
    are two or three regions, optimized C functions are used. Where
    there are more than three regions, the non-optimized python
    algorithm is used.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the co-segregation \
            of all possible combinations of windows within the different regions.
    """

    regions = prepare_regions(regions)

    if len(regions) == 2:
        coseg_func = cosegregation_2d
    elif len(regions) == 3:
        coseg_func = cosegregation_3d
    else:
        coseg_func = cosegregation_nd

    return coseg_func(*regions)


def get_cosesgregation(segregation_data, *location_strings):
    """Calculate co-segregation frequencies for a given genomic
    location or locations. Where only one location is given,
    co-segregation is calculated for that region against itself.

    :param segregation_data: Input :ref:`segregation table <segregation_table>`
    :param str location_strings: One or more :ref:`location strings <location_string>`
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the co-segregation \
            of all possible combinations of windows within the different regions.
    """

    regions = [segregation.region_from_location_string(
        segregation_data, l) for l in location_strings]

    return get_cosegregation_from_regions(*regions)


def get_linkage_from_regions(*regions):
    """Get the full linkage disequilibrium matrix for n regions, using optimized
    functions where available.

    This is a wrapper which determines the correct linkage
    function to call based on the number of regions. Where there
    are two regions, an optimized C function is used. Where
    there are three regions, linkage is calculated following the
    method given in
    `Hastings (1984) <http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1202243/>`_.
    This usage will raise a warning, since three-dimensional linkage is non-
    trivial and not necessarily obvious. Linkage matrices for more than three
    regions are not currently supported.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the linkage disequilibrium \
            of all possible combinations of windows within the different regions.
    """

    regions = prepare_regions(regions)

    if len(regions) == 2:
        linkage_func = linkage_2d
    elif len(regions) == 3:
        warnings.warn(
            '3D linkage is calculated following Hastings - Genetics (1984) 106:153-164. '
            'Please check and make sure this is what you want.')
        linkage_func = linkage_3d
    else:
        raise NotImplementedError(
            'No linkage function is defined for {} regions'.format(
                len(regions)))

    return linkage_func(*regions)


def get_linkage(segregation_data, *location_strings):
    """Calculate the linkage disequilibrium matrix for a given genomic
    location or locations. Where only one location is given,
    linkage is calculated for that region against itself.

    :param segregation_data: Input :ref:`segregation table <segregation_table>`
    :param str location_strings: One or more :ref:`location strings <location_string>`
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the linkage disequilibrium \
            of all possible combinations of windows within the different regions.
    """

    regions = [segregation.region_from_location_string(
        segregation_data, l) for l in location_strings]

    return get_linkage_from_regions(*regions)


def get_dprime_from_regions(*regions):
    """Get the full normalized linkage disequilibrium (D') matrix for n
    regions.

    This is a wrapper which determines the correct normalized linkage
    function to call based on the number of regions. Only two-dimensional
    normalized linkage matrices are currently supported. Where only one
    region is given, normalized linkage is calculated for that region against
    itself.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the normalized linkage \
            disequilibrium of all possible combinations of windows within the different regions.
    """

    regions = prepare_regions(regions)

    if len(regions) == 2:
        dprime_func = dprime_2d
    else:
        raise NotImplementedError(
            'There is currently no implementation of normalized linkage '
            'disequilibrium for more than 2 dimensions')

    return dprime_func(*regions)


def get_dprime(segregation_data, *location_strings):
    """Calculate the normalized linkage disequilibrium (D') matrix for a given
    genomic location or locations. Where only one location is given, normalized
    linkage is calculated for that region against itself.  Where two regions
    are given, linkage is calculated for region one against region two.

    :param segregation_data: Input :ref:`segregation table <segregation_table>`
    :param str location_strings: One or more :ref:`location strings <location_string>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the normalized linkage \
            disequilibrium of all possible combinations of windows within the different regions.
    """

    regions = [segregation.region_from_location_string(
        segregation_data, l) for l in location_strings]

    return get_dprime_from_regions(*regions)


def get_npmi_from_regions(*regions):
    """Get the full normalized pointwise mutual information (npmi) matrix for n
    regions.

    This is a wrapper which determines the correct npmi function to call based
    on the number of regions. Only two-dimensional matrices are currently
    supported. Where only one region is given, npmi is calculated for that
    region against itself.

    :param list regions: List of :ref:`regions <regions>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the \
            normalized pointwise mutual information (npmi) of all \
            possible combinations of windows within the different regions.
    """

    regions = prepare_regions(regions)

    if len(regions) == 2:
        npmi_func = npmi_2d
    else:
        raise NotImplementedError(
            'There is currently no implementation of npmi for more than '
            '2 dimensions')

    return npmi_func(*regions)


def get_npmi(segregation_data, *location_strings):
    """Calculate the normalized pointwise mutual information (npmi) matrix for a given
    genomic location or locations. Where only one location is given, npmi
    is calculated for that region against itself.  Where two regions
    are given, linkage is calculated for region one against region two.

    :param segregation_data: Input :ref:`segregation table <segregation_table>`
    :param str location_strings: One or more :ref:`location strings <location_string>`.
    :returns: :ref:`proximity matrix <proximity_matrices>` giving the \
            normalized pointwise mutual information (npmi) of all \
            possible combinations of windows within the different regions.
    """

    regions = [segregation.region_from_location_string(
        segregation_data, l) for l in location_strings]

    return get_npmi_from_regions(*regions)

MATRIX_TYPES = {
    'npmi': get_npmi_from_regions,
    'dprime': get_dprime_from_regions,
    'linkage': get_linkage_from_regions,
    'cosegregation': get_cosegregation_from_regions,
}


def get_regions_and_windows(segregation_data, location_strings):
    """Get the windows which fall into a given genomic location, and the
    segregation of those windows across samples.

    For each location string, find the rows of the :ref:`segregation table <segregation_table>`
    that overlap the given genomic location. Return the resulting
    :ref:`segregation table <segregation_table>` subsets (i.e. :ref:`regions`),
    and a list of the locations of the windows within each region.

    :param segregation_data: Input :ref:`segregation table <segregation_table>`
    :param list location_strings: One or more :ref:`location strings <location_string>`
    :returns: A list of :ref:`regions` and a list of tuples giving \
            window locations in the form (chromosome, start, stop).
    """

    if len(location_strings) == 1:
        location_strings = location_strings * 2

    regions = [
        segregation.region_from_location_string(
            segregation_data,
            location_str) for location_str in location_strings]

    windows = [np.array(list(region.index)) for region in regions]

    return regions, windows


def matrix_from_segregation_file(
        segregation_file, location_strings, matrix_type='dprime'):
    """Get the proximity matrix between the given genomic locations, and the
    locations of the genomic windows corresponding to each axis of the
    proximity matrix.

    Calculate the :ref:`proximity matrix <proximity_matrices>` between the
    windows in different genomic locations. If only one location string is
    given, calculate the matrix for that location against itself. For example,
    if location_strings is ['chr1'], the matrix will give the proximity between
    all windows on chromosome 1 (x-axis) against all windows on chromosome 1
    (y-axis).  Alternatively, if location_strings is ['chr1', 'chr2'], the
    matrix will give the proximity between windows on chromosome 1 (x-axis) and
    windows on chromosome 2 (y-axis).

    Three types of :ref:`proximity matrix <proximity_matrices>` are supported, 'cosegregation' (see
    :func:`get_cosegregation_from_regions`), 'linkage' (see
    :func:`get_linkage_from_regions`) and 'dprime' (see
    :func:`get_dprime_from_regions`).

    The function also returns a list of the windows that correspond to
    the axes of the :ref:`proximity matrix <proximity_matrices>`.

    :param segregation_file: Path to input :ref:`segregation table <segregation_table>`
    :param list location_strings: One or more :ref:`location strings <location_string>`
    :returns: A :ref:`proximity matrix <proximity_matrices>` for the given \
            genomic locations, and a list of tuples giving window locations \
            in the form (chromosome, start, stop).

    """

    segregation_data = segregation.open_segregation(segregation_file)

    regions, windows = get_regions_and_windows(
        segregation_data, location_strings)
    matrix_func = MATRIX_TYPES[matrix_type]
    contact_matrix = matrix_func(*regions)

    return contact_matrix, windows


def create_and_save_contact_matrix(segregation_file, location_strings,
                                   output_file, output_format,
                                   matrix_type='dprime'):
    """Calculate the proximity matrix for the given genomic locations and save it
    to disk.

    Calculate the :ref:`proximity matrix <proximity_matrices>` between the
    windows in different genomic locations. If only one location string is
    given, calculate the matrix for that location against itself. For example,
    if location_strings is ['chr1'], the matrix will give the proximity between
    all windows on chromosome 1 (x-axis) against all windows on chromosome 1
    (y-axis).  Alternatively, if location_strings is ['chr1', 'chr2'], the
    matrix will give the proximity between windows on chromosome 1 (x-axis) and
    windows on chromosome 2 (y-axis).

    Three types of :ref:`proximity matrix <proximity_matrices>` are supported, 'cosegregation' (see
    :func:`get_cosegregation_from_regions`), 'linkage' (see
    :func:`get_linkage_from_regions`) and 'dprime' (see
    :func:`get_dprime_from_regions`).

    Three output formats are supported: 'npz', 'txt' and 'csv'. 'txt' and
    'csv' outputs also support gzip compression ('txt.gz' and 'csv.gz'
    formats). See :ref:`matrix_formats` for more details.

    :param str segregation_file: Path to input :ref:`segregation table <segregation_table>`.
    :param list location_strings: One or more :ref:`location strings <location_string>`.
    :param str output_file: Path to use for saving output file.
    :param str output_format: Format to use when saving matrix. \
            (see :ref:`matrix_formats` for more details)
    :param str matrix_type: Type of :ref:`proximity matrix <proximity_matrices>`\
            to calculate.
    """

    print('starting calculation for {}'.format(' x '.join(location_strings)))
    start_time = time.perf_counter()

    contact_matrix, windows = matrix_from_segregation_file(
        segregation_file, location_strings, matrix_type)

    size_string = ' x '.join([str(s) for s in contact_matrix.shape])
    print('region size is: {}'.format(size_string), end=' ')
    print('Calculation took {0}s'.format(time.perf_counter() - start_time))
    print('Saving matrix to file {}'.format(output_file))

    output_func = matrix.OUTPUT_FORMATS[output_format]

    output_func(windows, contact_matrix, output_file)

    print('Done!')


def get_output_file(
        segregation_file,
        location_strings,
        matrix_type,
        output_format):
    """Automatically generate an output file path for a proximity matrix given
    the input segregation file.

    This function generates an output path in the format:

    <input_segregation_path>.<region1>_x_<region2>_<matrix_type>.<extension>

    :param str segregation_file: Path to input :ref:`segregation table <segregation_table>`.
    :param list location_strings: One or more :ref:`location strings <location_string>`.
    :param str matrix_type: Type of :ref:`proximity matrix <proximity_matrices>`\
            to calculate.
    :param str output_format: Format to use when saving matrix. \
            (see :ref:`matrix_formats` for more details)
    :returns: Path to save matrix file.

    >>> get_output_file('/path/to/segregation_file.table', ['chr1'], 'dprime', 'txt.gz')
    '/path/to/segregation_file.chr1_dprime.txt.gz'
    >>> get_output_file('/path/to/segregation_file.table', ['chr1', 'chr2'], 'linkage', 'npz')
    '/path/to/segregation_file.chr1_x_chr2_linkage.npz'
    >>> get_output_file('/path/to/segregation_file.table', ['chr1:10000000-20000000'],
    ...                 'cosegregation', 'csv')
    '/path/to/segregation_file.chr1_10.0Mb-20.0Mb_cosegregation.csv'

    """

    segregation_base = '.'.join(segregation_file.split('.')[:-1])
    locations = []
    for loc_string in location_strings:

        # If there is no colon, we are asked for a whole chromosome
        if not ':' in loc_string:
            locations.append(loc_string)
            continue

        chrom, start, stop = segregation.parse_location_string(loc_string)
        start, stop = format_genomic_distance(
            start), format_genomic_distance(stop)
        formatted_location = '{chrom}_{start}-{stop}'.format(chrom=chrom,
                                                             start=start,
                                                             stop=stop)
        locations.append(formatted_location)
    region_string = '_x_'.join(locations)

    return '{}.{}_{}.{}'.format(
        segregation_base,
        region_string,
        matrix_type,
        output_format)


def matrix_from_args(args):
    """Extract parameters from an argparse namespace object and pass them to
    create_and_save_contact_matrix.
    """

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
            args.segregation_file, args.regions,
            args.matrix_type, args.output_format)

    elif args.output_file == '-':
        args.output_file = sys.stdout

    create_and_save_contact_matrix(args.segregation_file, args.regions,
                                   args.output_file, args.output_format,
                                   args.matrix_type)


def matrix_from_doit(output_file, segregation_file, region):
    """Partial function that passes parameters from doit to
    create_and_save_contact_matrix
    """

    create_and_save_contact_matrix(segregation_file, region,
                                   output_file, 'txt.gz', 'dprime')
