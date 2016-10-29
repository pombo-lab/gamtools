"""
.. _segregation_table:

The segregation module
======================

The primary output of a GAM experiment is a segregation table. The segregation
module contains functions for reading segregation tables and for extracting
specific subsets of a segregation table.

Segregation Tables
------------------

The segregation table contains all the information about which genomic loci
are found in which nuclear profiles (NPs). In other words, the segregation
table gives the segregation of genomic loci across a collection of NPs. Each
column in a segregation table represents one NP. Each row represents a given
genomic locus, and each cell contains a 1 (indicating that the locus was
present in that NP) or a 0 (indicating that the locus was absent).  The
segregation table contains all the information about chromatin folding that was
captured during the experiment. For example, the nuclear proximity between any
two loci is estimated based on the number of times those two loci are
found in the same NP. Segregation tables are :ref:`pandas.Dataframe` objects.


.. _regions:

Regions
-------

A number of functions in the gamtools package take one or more `regions` as
an input. Regions are simply subsets of a segregation table, generally a subset
of the rows from a segregation table that span a smaller genomic region. Regions
are easily extracted from segregation tables by using the
`func:region_from_location_string` function.



"""
import numpy as np
import pandas as pd
import os


class InvalidChromError(Exception):
    """Exception to be raised when an invalid chromosome is specified"""
    pass


def open_segregation(path_or_buffer):

    return pd.read_csv(path_or_buffer,
                       index_col=[0, 1, 2],
                       delim_whitespace=True)


def index_from_interval(segregation_data, interval):

    chrom, start, stop = interval

    if not start < stop:
        raise ValueError(
            'Interval start {0} larger than interval end {1}'.format(
                *interval))

    window_in_region = np.logical_and(
        np.logical_and(
            segregation_data.index.get_level_values('stop') > start,
            segregation_data.index.get_level_values('start') < stop),
        segregation_data.index.get_level_values('chrom') == chrom)

    covered_windows = np.nonzero(window_in_region)[0]

    if not len(covered_windows):
        if not chrom in segregation_data.index.levels[0]:
            raise InvalidChromError(
                '{0} not found in the list of windows'.format(chrom))

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

        start, stop = (int(pos.translate(None, ",")) for pos in pos_fields)

    return chrom, start, stop


def index_from_location_string(segregation_data, loc_string):

    interval = parse_location_string(loc_string)

    return index_from_interval(segregation_data, interval)


def region_from_location_string(segregation_data, loc_string):

    ix_start, ix_stop = index_from_location_string(
        segregation_data, loc_string)

    return segregation_data.iloc[ix_start:ix_stop, ]


def detection_frequencies(region):
    return region.sum(axis=1).astype(float) / \
        region.count(axis=1).astype(float)


def map_sample_name_to_column(segregation_data):

    name_mapping = {}

    for c in segregation_data.columns:

        name_mapping[os.path.basename(c).split('.')[0]] = c

    return name_mapping


def sample_segregation_to_bed(input_segregation, sample, output_bed):

    data = open_segregation(input_segregation)

    subset = data[data[sample] > 0][sample]

    subset.to_csv(output_bed, header=False, index=True, sep='\t')


def is_autosome(chrom):
    if chrom[-7:] == '_random':
        return False
    try:
        int(chrom[3:])
    except ValueError:
        return False
    return True


def get_segregation_autosomes(segregation_table):
    return [c for c in segregation_table.index.get_level_values(0).unique()
            if is_autosome(c)]
