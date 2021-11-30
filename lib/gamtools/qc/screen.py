"""
====================
The qc.screen module
====================

The qc.screen module contains functions for parsing the output of
fastq_screen.

"""

import os

import pandas as pd

def screen_out_path(input_fastq):
    """
    Determine the path to the output file created by fastq_screen given
    an input fastq file.

    :param str input_fastq: Path to input fastq or fastq.gz
    :returns: Path to fastq_screen output file
    """

    path_parts = input_fastq.split('.')

    if path_parts[-1] == 'gz':
        path_parts.pop()

    if path_parts[-1] in ['txt', 'seq', 'fastq', 'fq']:
        path_parts.pop()

    path_parts[-1] += '_screen'
    path_parts.append('txt')

    return '.'.join(path_parts)

def get_sample_from_screen_path(fastq_screen_path):
    """
    Determine the name of the sample from the path to the fastq_screen output
    file.

    :param str fastq_screen_path: Path to fastq_screen output file
    :returns: Name of the sample
    """

    sample = os.path.basename(fastq_screen_path).split('.')[0]
    if sample[-7:] == '_screen':
        sample = sample[:-7]

    return sample

def is_fq_screen_header_row(fields):
    """
    Returns true if the input list represents a header row.
    """

    return bool((len(fields) == 0) or
                (fields[0][0] == '#') or
                (fields[0] == 'Library'))

def process_fastq_screen_line(line):
    """
    Split one line of a fastq_screen file on whitespace and parse the results.

    :param str line: One line of a fastq_screen output file
    :returns: Dictionary summarizing line results
    """

    fields = line.strip().split()

    if is_fq_screen_header_row(fields):
        row_results = {}

    elif fields[0] == '%Hit_no_libraries:':
        row_results = {'Unmapped': float(fields[1])}
    else:
        row_results = {
            fields[0] + '_single': int(fields[4]) + int(fields[8]),
            fields[0] + '_multiple': int(fields[6]) + int(fields[10]),
            'num_reads': int(fields[1]),
        }

    return row_results

def parse_fastq_screen_output(fastq_screen_output):
    """
    Parse the output of a single fastq_screen file.

    :param file fastq_screen_output: Open file object containing \
            fastq_screen results.
    :returns: Dictionary of percentage mapped reads for each \
            genome searched.
    """

    results = {}

    for line in fastq_screen_output:
        try:
            results.update(process_fastq_screen_line(line))
        except ValueError as error:
            raise ValueError('Malformed line: "{}"'.format(line)) from error


    total_reads = results['num_reads']
    del results['num_reads']

    for key, value in results.items():
        if key == 'Unmapped':
            continue
        results[key] = 100 * float(value) / total_reads

    organisms = []

    for key in results:

        try:
            organism, map_type = key.split('_')
        except ValueError:
            continue

        if map_type not in ['single', 'multiple']:
            continue

        if organism in organisms:
            continue

        organisms.append(organism)


    for organism in organisms:

        single_key, multi_key = organism + '_single', organism + '_multiple'

        results[organism] = results[single_key] + results[multi_key]

        del results[single_key]
        del results[multi_key]

    results = {key + '_screen_result': value for key, value in results.items()}

    return results


def get_contamination_stats(fastq_screen_output_files):
    """
    Process a list of fastq_screen files.

    :param list fastq_screen_output_files: List of paths to fastq_screen \
            output files.
    :returns: :class:`~pandas.DataFrame` containing percentage mapped reads \
            to each target genome for each NP.
    """

    sample_contamination = []

    for screen_file in fastq_screen_output_files:

        with open(screen_file) as screen_results:

            results = parse_fastq_screen_output(screen_results)

            results['Sample'] = get_sample_from_screen_path(screen_file)

            sample_contamination.append(results)

    contam_df = pd.DataFrame(sample_contamination)

    columns = ['Sample'] + [col for col in contam_df.columns if col != 'Sample'] #pylint: disable=not-an-iterable

    return contam_df[columns]

def write_contamination_stats(input_files, output_file):
    """
    Given a list of fastq_screen output files, parse the files to a
    dataframe and save the dataframe to output file.

    :param list input_files: List of paths to fastq_screen \
            results files.
    :param str output_file: Path to save output file.
    """

    contam_df = get_contamination_stats(input_files)

    contam_df.to_csv(output_file, sep='\t', index=False)

def contamination_from_doit(dependencies, targets):
    """Wrapper function to call write_contamination_stats from argparse"""

    assert len(targets) == 1
    write_contamination_stats(dependencies, targets[0])
