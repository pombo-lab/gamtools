"""This module holds utility functions that can be called
by other components but don't have their own place.
"""


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

