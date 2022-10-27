"""
================
The utils module
================

The utils module contains various utility functions that don't have
their own place.

"""

import sys
import os

import pandas as pd


def get_example(example_name):
    """Get the absolute path to a file in the gamtools 'examples' folder

    :param str example_name: Name of the file.
    :returns: Absolute path to the file
    """

    return os.path.join(os.path.dirname(__file__),
                        'data',
                        'examples',
                        example_name)


def format_genomic_distance(distance, precision=1):
    """
    Turn an integer genomic distance into a pretty string.

    :param int distance: Genomic distance in basepairs.
    :param int precision: Number of significant figures to display \
        after the decimal point.
    """

    formatting_string = '{{0:.{0}f}}'.format(precision)

    if distance < 1000:
        return '{0:d}bp'.format(int(distance))

    if distance < 1000000:
        fmt_string = formatting_string + 'kb'
        return fmt_string.format(float(distance) / 1000)

    fmt_string = formatting_string + 'Mb'
    return fmt_string.format(float(distance) / 1000000)


def empty_bedgraph(chrom_sizes, output_bedgraph):
    """
    Create an empty bedgraph file.

    :param str chrom_sizes: Path to file containing a list of chromosome names and sizes
    :param str output_bedgraph: Path to save bedgraph.
    """

    bed_df = pd.read_csv(chrom_sizes, sep='\t', names=['chrom', 'stop'])

    bed_df['start'] = 1
    bed_df['score'] = 0

    bed_df[['chrom', 'start', 'stop', 'score']].to_csv(
        output_bedgraph, sep='\t', index=False, header=False)


def empty_bedgraph_from_cmdline():
    """
    Wrapper function to create empty bedgraphs from the command line.
    """

    empty_bedgraph(sys.argv[1], sys.argv[2])


class DelayedImportError(): #pylint: disable=too-few-public-methods
    """Class that returns an ImportError if any method or attribute is accessed.

    Useful for delaying the ImportError until an optional dependency is actually used"""

    def __init__(self, message):
        """Instantiate the object with an error message and a list
        of packages that need to be installed"""

        self.message = message

    def __getattr__(self, name):
        """Raise ImportError if any method or attribute is called"""
        raise ImportError(self.message)

    def __call__(self, *args, **kwargs):
        """Raise ImportError if the object itself is called"""
        raise ImportError(self.message)
