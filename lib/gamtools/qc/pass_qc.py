"""
====================
The qc.passqc module
====================

The qc.passqc module contains functions for determining which NPs pass
a set of quality control conditions.

"""

import pandas as pd


def get_reference(condition_ref_string, stats_df):
    """
    If condition_ref_string matches a column in stats_df, return that
    column. If it doesn't match a column, try to convert it to a float.

    :param str condition_ref_string: left or right part of a condition \
            from a QC conditions file
    :param stats_df: QC statistics table
    :type stats_df: :class:`~pandas.DataFrame`
    :returns: Either a column of stats_df, a float or a string.
    """

    if condition_ref_string in stats_df.columns:
        return stats_df.loc[:, condition_ref_string]

    try:
        return float(condition_ref_string)
    except ValueError:
        return condition_ref_string


def get_references(left_str, right_str, stats_df):
    """
    Parse the left and right references in one line of a QC conditions file.

    e.g. "left > right" or "left less_than 5.0"

    :param str left_str: Left part of a condition \
            from a QC conditions file
    :param str right_str: Right part of a condition \
            from a QC conditions file
    :param stats_df: QC statistics table
    :type stats_df: :class:`~pandas.DataFrame`
    :returns: left reference, right reference (at least one must be a column of \
            stats_df. The other can be another column, a float or a string)
    """

    left = get_reference(left_str, stats_df)
    right = get_reference(right_str, stats_df)

    if not (isinstance(left, pd.Series) or isinstance(right, pd.Series)):
        raise QcParamError(
            """Neither {} nor {} matches a column in the stats file.
            Current options are: \n\t'{}'""".format(
                left_str, right_str,
                "'\n\t'".join(stats_df.columns)))

    return left, right


def comparison_from_operator(operator, left, right):
    """
    Perform a comparison between left and right values in a QC
    conditions file.

    :param str operator: Comparison to carry out.
    :param left: Either a series of QC values or a value to compare QC values to.
    :param right: Either a series of QC values or a value to compare QC values to.
    :returns: :class:`~pandas.Series` indicating which samples pass the condition.
    """

    if operator in ['=', '==', 'eq', 'equals']:
        comparison = (left == right)
    elif operator in ['>', 'gt', 'greater_than']:
        comparison = (left > right)
    elif operator in ['>=', 'gte', 'greater_than_or_equal_to']:
        comparison = (left >= right)
    elif operator in ['<', 'lt', 'less_than']:
        comparison = (left < right)
    elif operator in ['<=', 'lte', 'less_than_or_equal_to']:
        comparison = (left <= right)
    elif operator in ['!=', 'neq', 'not_equal_to']:
        comparison = (left != right)
    else:
        raise QcParamError(
            'Operator {} not recognized'.format(operator))

    return comparison


class QcParamError(Exception):
    """
    Exception to be raised if the QC Parameters file is malformed.
    """


def do_comparison(left_str, operator, right_str, stats_df):
    """
    Given the condition in a QC Parameters file, and a QC statistics table,
    calculate which samples pass the condition.

    :param str left_str: Left part of the condition (e.g. "mapped_reads")
    :param str operator: How to compare left and right (e.g. greater_than)
    :param str right_str: Right part of the condition (e.g. "150000")
    :param stats_df: QC statistics table
    :type stats_df: :class:`~pandas.DataFrame`
    :returns: :class:`~pandas.Series` indicating which samples pass the condition.
    """

    left, right = get_references(left_str, right_str, stats_df)

    comparison = comparison_from_operator(operator, left, right)

    if not isinstance(comparison, pd.Series):
        raise QcParamError('Comparison did not return a series object')

    return comparison


def parse_conditions_file(conditions_file, stats_df):
    """
    Iterate through lines of a conditions file and perform
    the indicated comparison for each line.

    :param file conditions_file: Open file buffer containing conditions.
    :param stats_df: QC statistics table
    :type stats_df: :class:`~pandas.DataFrame`
    :returns: List of :class:`~pandas.Series` each indicating whether each \
            NP passed each quality check.
    """

    conditions = []

    for line_no, line in enumerate(conditions_file, 1):

        fields = line.split()

        if (line[0] == '#') or (not fields):
            continue

        left_str, operator, right_str = fields

        try:
            this_comparison = do_comparison(left_str, operator, right_str, stats_df)
        except QcParamError as err_msg:
            raise QcParamError(
                'Error in QC conditions file line {}: {}'.format(
                    line_no, err_msg)) from err_msg

        conditions.append(this_comparison)

    return conditions


def samples_passing_qc(conditions_file_path, stats_df_path):
    """
    Read a QC conditions file and a QC stats file, calculating which
    samples pass all specified QC checks.

    :param str conditions_file_path: Path to QC conditions file.
    :param stats_df_path: Path to QC statistics table
    :returns: :class:`~pandas.Series` of samples that pass all specified \
            QC checks.
    """

    stats_df = pd.read_csv(stats_df_path, delim_whitespace=True)

    with open(conditions_file_path, 'r') as conditions_file:
        conditions = parse_conditions_file(conditions_file, stats_df)

    sample_passes_qc = pd.concat(conditions, axis=1).all(axis=1)

    return stats_df.loc[sample_passes_qc, 'Sample']

def create_passqc_file(conditions_file_path, stats_df_path, passqc_path):
    """
    Read a QC conditions file and a QC stats file, calculate which
    samples pass all specified QC checks and save that list of samples
    to a file.

    :param str conditions_file_path: Path to QC conditions file.
    :param stats_df_path: Path to QC statistics table
    :param passqc_path: Path to save IDs of samples that pass all QC \
            checks.
    """


    samples_passing = samples_passing_qc(conditions_file_path, stats_df_path)
    samples_passing = samples_passing.sort_values()
    samples_passing.to_csv(passqc_path, sep='\t', index=False, header=False)

def samples_passing_qc_from_doit(targets, dependencies):
    """
    Wrapper function to call create_passqc_file from argparse.
    """

    for file_path in dependencies:

        file_ext = file_path.split('.')[-1]
        if file_ext == 'cfg':
            conditions_file_path = file_path
        elif file_ext == 'txt':
            stats_df_path = file_path

    create_passqc_file(conditions_file_path, stats_df_path, list(targets)[0])
