"""This module holds utility functions that can be called
by other components but don't have their own place.
"""

import pandas as pd
import sys


def format_genomic_distance(distance, precision=1):
    """Turn an integer genomic distance into a pretty string.
    :param int distance: Genomic distance in basepairs.
    :param int precision: Number of significant figures to display
        after the decimal point.
    """

    formatting_string = '{{0:.{0}f}}'.format(precision)

    if distance < 1000:
        return '{0:d}bp'.format(int(distance))
    elif distance < 1000000:
        fmt_string = formatting_string + 'kb'
        return fmt_string.format(float(distance) / 1000)
    else:
        fmt_string = formatting_string + 'Mb'
        return fmt_string.format(float(distance) / 1000000)

def empty_bedgraph(chrom_sizes, output_bedgraph):

    bed_df = pd.read_csv(chrom_sizes, sep='\t', names=['chrom', 'stop'])

    bed_df['start'] = 1
    bed_df['score'] = 0

    bed_df[['chrom', 'start', 'stop', 'score']].to_csv(output_bedgraph, sep='\t', index=False, header=False)

def empty_bedgraph_from_cmdline():

    empty_bedgraph(sys.argv[1], sys.argv[2])
