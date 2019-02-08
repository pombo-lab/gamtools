"""

The segregation module
======================

The primary output of a GAM experiment is a segregation table. The segregation
module contains functions for reading segregation tables and for extracting
specific subsets of a segregation table.

.. _NPs:

Nuclear profiles
----------------

In a GAM experiment, thin sections are taken through a cell nucleus and
the DNA content is sequenced. Each slice of a single nucleus is called
a nuclear profile or "NP".

.. _segregation_table:

Segregation tables
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
found in the same NP. Segregation tables are :class:`pandas.DataFrame` objects.


.. _regions:

Regions
-------

A number of functions in the gamtools package take one or more `regions` as
an input. Regions are simply subsets of a segregation table, generally a subset
of the rows from a segregation table that span a smaller genomic region. Regions
are easily extracted from segregation tables by using the
`func:region_from_location_string` function.

.. _location_string:

Location strings
----------------

UCSC_ style location strings are convenient ways to specify a given genomic
region. They are given in the format "chromosome:start-end", e.g.
"chr2:100000-200000" or "chrX:15001-15100". The names of different chromosomes
are typically determined by the genome annotation you are using.

.. _UCSC: https://genome.ucsc.edu/


"""
import os

import numpy as np
import pandas as pd


class InvalidChromError(Exception):
    """Exception to be raised when an invalid chromosome is specified"""


def open_segregation(path_or_buffer):
    """
    Open a segregation table from a file.

    :param path_or_buffer: Path to imput segregation table, or open python file object.
    :returns: :ref:`segregation table <segregation_table>`
    """

    return pd.read_csv(path_or_buffer,
                       index_col=[0, 1, 2],
                       delim_whitespace=True)


def index_from_interval(segregation_table, interval):
    """
    Find the start and end rows of the sub-region of a segregation table.
    The region is given as a tuple in the form (chromosome, start_coordinate, end_coordinate).

    :param segregation_table: Input segregation table.
    :type segregation_table: :ref:`segregation table <segregation_table>`
    :param tuple interval: Region of segregation table to search for.
    :returns: First row overlapping region, last row overlapping region
    """

    chrom, start, stop = interval

    if not start < stop:
        raise ValueError(
            'Interval start {0} larger than interval end {1}'.format(
                *interval))

    window_in_region = np.logical_and( #pylint: disable=assignment-from-no-return
        np.logical_and(
            segregation_table.index.get_level_values('stop') > start,
            segregation_table.index.get_level_values('start') < stop),
        segregation_table.index.get_level_values('chrom') == chrom)

    covered_windows = np.nonzero(window_in_region)[0]

    if not covered_windows.size:
        if not chrom in segregation_table.index.levels[0]:
            raise InvalidChromError(
                '{0} not found in the list of windows'.format(chrom))

    start_index = covered_windows[0]
    stop_index = covered_windows[-1] + 1

    return start_index, stop_index


def parse_location_string(loc_string):
    """
    Parse a UCSC format location string (e.g. "chr2:1000-1100") and return
    an interval tuple in the format ('chr2', 1000, 1100).

    :param loc_string: Input location string
    :type loc_string: :ref:`location string <location_string>`
    :returns: (chromosome name, start coordinate, stop coordinate)
    """

    chrom_fields = loc_string.split(':')

    chrom = chrom_fields[0]

    if len(chrom_fields) == 1:

        start, stop = 0, np.Inf

    else:

        pos_fields = chrom_fields[1].split('-')

        start, stop = (int(pos.replace(",", "")) for pos in pos_fields)

    return chrom, start, stop


def index_from_location_string(segregation_table, loc_string):
    """
    Find the start and end rows of the sub-region of a segregation table.
    The region is given in the form of a UCSC format location string
    (e.g. "chr2:1000-1100").

    :param segregation_table: Input segregation table.
    :type segregation_table: :ref:`segregation table <segregation_table>`
    :param loc_string: Input location string
    :type loc_string: :ref:`location string <location_string>`
    :returns: First row overlapping region, last row overlapping region
    """

    interval = parse_location_string(loc_string)

    return index_from_interval(segregation_table, interval)


def region_from_location_string(segregation_table, loc_string):
    """
    Extract the rows of a segregation table that overlap a given genomic region.
    The region is given in the form of a UCSC format location string
    (e.g. "chr2:1000-1100").

    :param segregation_table: Input segregation table.
    :type segregation_table: :ref:`segregation table <segregation_table>`
    :param loc_string: Input location string
    :type loc_string: :ref:`location string <location_string>`
    :returns: Subset of the input segregation table that overlaps the region.
    """

    ix_start, ix_stop = index_from_location_string(
        segregation_table, loc_string)

    return segregation_table.iloc[ix_start:ix_stop, ]


def detection_frequencies(segregation_table):
    """
    Calculate the detection frequency of each window in a segregation table
    across NPs. Detection frequency is the number of NPs a window was detected
    in divided by the total number of NPs.

    :param segregation_table: Input segregation table.
    :type segregation_table: :ref:`segregation table <segregation_table>`
    :returns: Pandas Series object containin detection frequency of each window
    """

    return segregation_table.sum(axis=1).astype(float) / \
        segregation_table.count(axis=1).astype(float)


def map_sample_name_to_column(segregation_table):
    """
    Return a dictionary mapping NP sample names to column names in a
    segregation table.

    :param segregation_table: Input segregation table.
    :type segregation_table: :ref:`segregation table <segregation_table>`
    :returns: Dictionary of sample_name:column_name mappings
    """

    name_mapping = {}

    for column in segregation_table.columns:

        name_mapping[os.path.basename(column).split('.')[0]] = column

    return name_mapping


def sample_segregation_to_bed(input_segregation, sample, output_bed):
    """
    Create a bed file containing the positive windows for one sample in a
    segregation table.

    :param str input_segregation: Path to imput segregation table, or open python file object.
    :param str sample: Column name for which to extract positive windows
    :param str output_bed: Path to save bed file.
    """

    data = open_segregation(input_segregation)

    subset = data[data[sample] > 0][sample]

    subset.reset_index()[['chrom', 'start', 'stop']].to_csv(
        output_bed, header=False, index=False, sep='\t')


def is_autosome(chrom):
    """
    Returns true if the given chromosome label likely refers to an autosome.

    :param str chrom: Chromosome label
    """

    if chrom[-7:] == '_random':
        return False
    try:
        int(chrom[3:])
    except ValueError:
        return False
    return True


def get_segregation_autosomes(segregation_table):
    """
    Return the autosomal chromosomes listed in a segregation table.

    :param segregation_table: Input segregation table.
    :type segregation_table: :ref:`segregation table <segregation_table>`
    :returns: List of autosomal chromosomes in the segregation table.
    """


    return [c for c in segregation_table.index.get_level_values(0).unique()
            if is_autosome(c)]
