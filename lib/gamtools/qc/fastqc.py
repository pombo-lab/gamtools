"""
=================
The qc.fastqc module
=================

The qc.fastqc module contains functions for parsing fastqc output files.
Code in this module was shamelessly stolen from
https://code.google.com/p/bioinformatics-misc/source/browse/trunk/fastqc_to_pgtable.py?spec=svn93&r=93

"""


import os

import numpy as np
import pandas as pd


def fastqc_data_file(input_fastq):
    """Given an input fastq file, return the name fastqc will use
    for it's output file.

    :param str input_fastq: Path to a fastq file.
    :returns: Path to fastqc output files.
    """

    base_folder = input_fastq.split('.')[0]
    #base_folder = '.'.join(input_fastq.split('.')[:-1])

    fastqc_folder = base_folder + '_fastqc'

    return os.path.join(fastqc_folder, 'fastqc_data.txt')


def parse_module(fastqc_module):
    """
    Parse a fastqc module from the table format to a line format (list).
    Input is list containing the module. One list-item per line. E.g.:

    fastqc_module= [
        '>>Per base sequence quality    pass',
        '#Base  Mean    Median  Lower Quartile  Upper Quartile  10th Percentile 90th Percentile',
        '1  36.34   38.0    35.0    39.0    31.0    40.0',
        '2  35.64   38.0    35.0    39.0    28.0    40.0',
        '3  35.50   38.0    34.0    39.0    28.0    40.0',
        ...
        ]

    Return a list like this where each sublist after 1st is a column:
    ['pass', ['1', '2', '3', ...], ['36,34', '35.64', '35.50', ...], ['40.0', '40.0', '40.0', ...]]
    """
    row_list = []
    module_header = fastqc_module[0]
    module_name = module_header.split('\t')[0]
    # Line with module name >>Per base ...    pass/warn/fail
    row_list.append(module_header.split('\t')[1])

    # Handle odd cases:
    # Where no table is returned:
    if len(fastqc_module) == 1 and module_name == '>>Overrepresented sequences':
        return row_list + [[]] * 4
    if len(fastqc_module) == 1 and module_name == '>>Kmer Content':
        return row_list + [[]] * 5
    # Table is not the secod row:
    if module_name == '>>Sequence Duplication Levels':
        tot_dupl = fastqc_module[1].split('\t')[1]
        row_list.append(tot_dupl)
        del fastqc_module[1]

    # Conevrt table to list of lists:
    tbl = []
    for line in fastqc_module[2:]:
        tbl.append(line.split('\t'))
    # Put each column in a list:
    nrows = len(tbl)
    ncols = len(tbl[0])
    for i in range(0, ncols):
        col = []
        for j in range(0, nrows):
            col.append(tbl[j][i])
        row_list.append(col)
    return row_list


def is_mono_repeat(kmer):
    """
    Return true if kmer represents a mono-nucleotide repeat (e.g. AAAA).
    """

    return bool(len(set(kmer)) == 1)


def is_di_repeat(kmer):
    """
    Return true if kmer represents a di-nucleotide repeat (e.g. ATATAT).
    """

    if not len(set(kmer)) == 2:
        return False

    return bool(len(set(kmer[::2])) == 1 and len(set(kmer[1::2])) == 1)


def get_kmer_summary(module):
    """
    Calculate the total number of kmers that are either mono- or di- nucleotide
    repeats.
    """

    kmer_data = parse_module(module)
    kmers = kmer_data[1]
    counts = list(map(float, kmer_data[3]))

    summary_data = {'dinucleotide_repeats': 0,
                    'mononucleotide_repeats': 0}

    for kmer, count in zip(kmers, counts):
        if is_mono_repeat(kmer):
            summary_data['mononucleotide_repeats'] += count
        elif is_di_repeat(kmer):
            summary_data['dinucleotide_repeats'] += count

    return summary_data


def get_avg_qual(module):
    """
    Get the average per-base-pair sequencing quality score.
    """

    qual_data = parse_module(module)
    qualities, counts = np.array(
        list(map(int, qual_data[1]))), np.array(list(map(float, qual_data[2])))
    avg = (qualities * counts).sum() / counts.sum()
    return {'avg_quality': avg}


def get_sample(filename):
    """
    Get the name of the input sample given the fastqc output file path.
    """

    return os.path.basename(os.path.dirname(filename))[:-7]


def process_file(filename):
    """
    Process a fastqc output file and calculate some summary statistics.
    """

    fq_lines = open(filename).readlines()
    fq_lines = [x.strip() for x in fq_lines]

    fastqc_dict = {}

    # Get start and end position of all modules:
    mod_start = []
    mod_end = []
    for i, line in enumerate(fq_lines):

        if line == '>>END_MODULE':
            mod_end.append(i)
        elif line.startswith('>>'):
            mod_start.append(i)
        else:
            pass

    # Start processing modules. First one (Basic statitics) is apart:
    for start, end in zip(mod_start[1:], mod_end[1:]):
        module = fq_lines[start:end]
        module_name = module[0].split('\t')[0]
        if module_name == '>>Kmer Content':

            fastqc_dict.update(get_kmer_summary(module))

        if module_name == '>>Per sequence quality scores':

            fastqc_dict.update(get_avg_qual(module))

    fastqc_dict['Sample'] = get_sample(filename)

    return fastqc_dict


def get_quality_stats(input_fastqc_files):
    """
    Iterate over a list of fastqc output files and generate a dataframe
    containing summary statistics for each file.
    """

    sample_qualities = []
    for filename in input_fastqc_files:
        sample_qualities.append(process_file(filename))

    return pd.DataFrame(sample_qualities)


def write_quality_stats(input_files, output_file):
    """
    Iterate over a list of fastqc output files and generate a dataframe
    containing summary statistics for each file, then write the result
    to disk.
    """

    quality_df = get_quality_stats(input_files)

    quality_df.to_csv(output_file, sep='\t', index=False)


def quality_qc_from_doit(dependencies, targets):
    """
    Wrapper function to call write_quality_stats from argparse.
    """

    assert len(targets) == 1
    write_quality_stats(dependencies, targets[0])
