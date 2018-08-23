import os
import errno

import numpy as np
import pandas as pd

from . import slice_wrapper, segregation

def mkdir(path):
    try:
        print("Making directory: {}".format(path))
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def get_slice_output_dirs(output_dir):

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
    print('Creating chr_indices file: {}'.format(chr_indices_path))
    bin_indices = segregation_table.reset_index().loc[:, ['chrom', 'start', 'stop']].reset_index()
    chr_indices = pd.concat([bin_indices.groupby('chrom').index.min().sort_values(),
                             bin_indices.groupby('chrom').index.max().sort_values()],
                            axis=1)
    chr_indices.to_csv(chr_indices_path, header=False, index=False, sep='\t')

def create_chr_names_file(segregation_table, chr_names_path):
    print('Creating chr_names file: {}'.format(chr_names_path))
    np.savetxt(chr_names_path,
               segregation_table.index.get_level_values(0).unique().values,
               fmt='%s')

def create_matrix_file(segregation_table, matrix_path):
    print('Creating matrix file: {}'.format(matrix_path))
    segregation_table.to_csv(matrix_path, header=False, index=False, sep='\t')

def segregation_info(segregation_table, haploid_chroms):
    m = segregation_table.shape[1]
    print('{} tubes in dataset'.format(m))

    segregation_windows = segregation_table.reset_index()[['chrom', 'start', 'stop']]
    segregation_windows['size'] = segregation_windows.stop - segregation_windows.start
    L = segregation_windows.loc[
            np.logical_not(segregation_windows.chrom.isin(haploid_chroms)),
            'size'].sum() * 2
    print('Diploid genome size is {}'.format(L))

    bin_sizes = segregation_windows['size'].value_counts()
    frequent_sizes = [size for (size, freq)
                      in (bin_sizes / bin_sizes.sum()).iteritems()
                      if freq > 0.99]
    if not frequent_sizes:
        raise Exception('SLICE requires windows of even sizes')
    b = frequent_sizes[0]
    print('Bin size is {}'.format(b))

    return m, L, b

def split_segregation_table(segregation_table, matrix_path, chr_names_path, chr_indices_path):

    create_chr_indices_file(segregation_table, chr_indices_path)
    create_chr_names_file(segregation_table, chr_names_path)
    create_matrix_file(segregation_table, matrix_path)

def process_segregation_table(segregation_file_path, matrix_path,
                              chr_names_path, chr_indices_path,
                              haploid_chroms):

    segregation_table = segregation.open_segregation(segregation_file_path)

    split_segregation_table(segregation_table, matrix_path, chr_names_path, chr_indices_path)

    return segregation_info(segregation_table, haploid_chroms)

def run_slice(segregation_file_path, output_dir, haploid_chroms, h, R):
    print('Running slice...')

    pi_out_path, threshold_out_path, chr_names_path, chr_indices_path, matrix_path = get_slice_output_dirs(output_dir)

    m, L, b = process_segregation_table(segregation_file_path, matrix_path,
                                     chr_names_path, chr_indices_path,
                                     haploid_chroms)

    slice_wrapper.slice(matrix_path,
          pi_out_path, threshold_out_path,
          chr_names_path, chr_indices_path,
          m, L, b, h, R)

def run_slice_from_args(args):
    run_slice(args.segregation_file_path,
              args.output_dir,
              args.haploid_chroms,
              args.slice_thickness,
              args.nuclear_radius)
