"""

The gam_slice module
======================

The SLICE module gives the user access to wrapper functions that allow the SLICE
library, which is written in C++ to be accessed from Python.

These functions allow the user to calculate "Probability of interaction" or PI
matrices from GAM segregation data.

"""

# Some variable names are inherited from the SLICE library
# pylint: disable=invalid-name
import os
import errno

import numpy as np
import pandas as pd

from . import slice_wrapper, segregation


def mkdir(path):
    """Helper function for easily making slice output directories"""

    try:
        print("Making directory: {}".format(path))
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def get_slice_output_dirs(output_dir):
    """
    Helper function that ensures all neccessary SLICE output directories exist
    before SLICE is run.
    """

    pi_out_dir = os.path.join(output_dir, 'pi')
    mkdir(pi_out_dir)
    pi_out_path = os.path.join(pi_out_dir, 'pi_pairs')

    threshold_out_dir = os.path.join(output_dir, 'pi_thresholds')
    mkdir(threshold_out_dir)
    threshold_out_path = os.path.join(threshold_out_dir, 'pi_thresholds')

    chr_details_dir = os.path.join(output_dir, 'chr_details')
    mkdir(chr_details_dir)
    chr_names_path = os.path.join(chr_details_dir, 'chr_names.out')
    chr_indices_path = os.path.join(chr_details_dir, 'chr_indices.out')

    matrix_path = os.path.join(output_dir, 'matrix.out')

    return pi_out_path, threshold_out_path, chr_names_path, chr_indices_path, matrix_path


def create_chr_indices_file(segregation_table, chr_indices_path):
    """
    Function that takes a segregation table and outputs a file containing the start
    and stop indexes of each chromosome
    """

    print('Creating chr_indices file: {}'.format(chr_indices_path))
    bin_indices = segregation_table.reset_index().loc[:, ['chrom', 'start', 'stop']].reset_index()
    chr_indices = pd.concat([bin_indices.groupby('chrom').index.min().sort_values(),
                             bin_indices.groupby('chrom').index.max().sort_values()],
                            axis=1)
    chr_indices.to_csv(chr_indices_path, header=False, index=False, sep='\t')


def create_chr_names_file(segregation_table, chr_names_path):
    """
    Given a segregation table, output a text file containing the chromosome names
    """

    print('Creating chr_names file: {}'.format(chr_names_path))
    np.savetxt(chr_names_path,
               segregation_table.index.get_level_values(0).unique().values,
               fmt='%s')


def create_matrix_file(segregation_table, matrix_path):
    """
    Save a segregation table in the correct input format for slice
    """

    print('Creating matrix file: {}'.format(matrix_path))
    segregation_table.to_csv(matrix_path, header=False, index=False, sep='\t')


def segregation_info(segregation_table, skip_chroms):
    """
    Given a segregation table, retrieve the number of tubes,
    the diploid genome size and the bin size.
    """

    m = segregation_table.shape[1]
    print('{} tubes in dataset'.format(m))

    segregation_windows = segregation_table.reset_index()[['chrom', 'start', 'stop']]
    segregation_windows['size'] = segregation_windows.stop - segregation_windows.start
    L = segregation_windows.loc[
        np.logical_not(segregation_windows.chrom.isin(skip_chroms)),
        'size'].sum() * 2

    bin_sizes = segregation_windows['size'].value_counts()
    frequent_sizes = [size for (size, freq)
                      in (bin_sizes / bin_sizes.sum()).iteritems()
                      if freq > 0.95]
    if not frequent_sizes:
        raise Exception('SLICE requires windows of even sizes')
    b = frequent_sizes[0]
    print('Bin size is {}'.format(b))

    return m, L, b


def split_segregation_table(segregation_table, matrix_path, chr_names_path, chr_indices_path):
    """
    Split the information in a segregation table into three different input files
    as required by the SLICE algorithm.
    """

    create_chr_indices_file(segregation_table, chr_indices_path)
    create_chr_names_file(segregation_table, chr_names_path)
    create_matrix_file(segregation_table, matrix_path)


def process_segregation_table(segregation_file_path, matrix_path,
                              chr_names_path, chr_indices_path,
                              skip_chroms):
    """
    Given a path to a segregation table, split the table into three input
    files required by SLICE, then return information about the table.
    """

    segregation_table = segregation.open_segregation(segregation_file_path)

    split_segregation_table(segregation_table, matrix_path, chr_names_path, chr_indices_path)

    return segregation_info(segregation_table, skip_chroms)


def run_slice(segregation_file_path, output_dir, skip_chroms, h, R, genome_size, n_p): #pylint: disable=too-many-arguments
    """
    Wrapper function that takes a path to a segregation file, opens and
    processes the segregation table and runs the underlying SLICE C++
    library.
    """

    print('Running slice...')

    (pi_out_path,
     threshold_out_path,
     chr_names_path,
     chr_indices_path,
     matrix_path) = get_slice_output_dirs(output_dir)

    m, L, b = process_segregation_table(segregation_file_path,
                                        matrix_path,
                                        chr_names_path,
                                        chr_indices_path,
                                        skip_chroms)

    if genome_size is not None:
        L = genome_size

    print('Slice thickness is {}'.format(h))
    print('Nuclear radius is {}'.format(R))
    print('Diploid genome size is {}'.format(L))
    print('{} NPs per tube'.format(n_p))

    slice_wrapper.slice(matrix_path, pi_out_path,
                        threshold_out_path,
                        chr_names_path, chr_indices_path,
                        m, L, b, h, R, n_p)


def run_slice_from_args(args):
    """Allow the run_slice function to be called from doit."""

    run_slice(args.segregation_file_path,
              args.output_dir,
              args.skip_chroms,
              args.slice_thickness,
              args.nuclear_radius,
              args.genome_size,
              args.nps_per_tube)
