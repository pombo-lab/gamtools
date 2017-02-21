"""
=========================
The select_samples module
=========================

The select_samples module contains functions for extracting samples
from :ref:`segregation_tables`.

"""

from . import segregation


def select_samples(segregation_path, sample_names, output_file, drop=False):
    """Create a file containing a subset of columns from a segregation table

    Given the path to a segregation table and a list of columns in that table,
    save a new version of the table containing only those columns to the path
    specified by output_file. If drop is True, save a new version of the table
    that excludes the listed columns, rather than including them

    :param str segregation_path: Path to input segregation table
    :param str sample_names: List of samples to subset on
    :param str segregation_path: Path to save output segregation table
    :param bool drop: Whether to exclude (drop) samples instead of including them.
    """

    segregation_data = segregation.open_segregation(segregation_path)

    column_mapping = segregation.map_sample_name_to_column(segregation_data)
    column_names = [column_mapping[name] for name in sample_names]

    if drop:
        subset = segregation_data.drop(column_names, axis=1)
    else:
        subset = segregation_data[column_names]

    subset.to_csv(output_file, index=True, sep='\t')


def select_samples_from_file(
        segregation_path,
        sample_names_path,
        output_file,
        drop=False):
    """Create a file containing a subset of columns from a segregation table

    Given the path to a segregation table and a path to a second file containing
    samples to subset on, save a new version of the table to the path
    specified by output_file. If drop is True, save a new version of the table
    that excludes the listed columns, rather than including them

    :param str segregation_path: Path to input segregation table
    :param str sample_names_path: Path to the list of samples to subset on
    :param str segregation_path: Path to save output segregation table
    :param bool drop: Whether to exclude (drop) samples instead of including them.
    """

    names = []

    with open(sample_names_path, 'r') as sample_names_file:
        for line in sample_names_file:
            names.append(line.strip())

    select_samples(segregation_path, names, output_file, drop)


def select_samples_from_args(args):
    """Wrapper function to call select_samples_from_file from argparse"""

    select_samples(
        args.segregation_file,
        args.sample_names,
        args.output_file,
        args.drop_samples)

def select_samples_from_doit(dependencies, targets):
    """Wrapper function to call select_samples_from_file from a doit task"""

    for dep_file in dependencies:
        dep_ext = dep_file.split('.')[-1]

        if dep_ext == 'table':
            segregation_path = dep_file
        elif dep_ext == 'txt':
            sample_names_path = dep_file

    assert len(targets) == 1

    select_samples_from_file(
        segregation_path,
        sample_names_path,
        list(targets)[0])
