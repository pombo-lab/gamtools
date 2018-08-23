"""
===============
The main module
===============

The `main` module creates the relatively complex argument parser used by the
GAMtools command line application, and the :func:`main` function which handles
dispatching command line arguments to the correct functions.

"""

# more or less everything in this module is in the global scope, so pylint
# thinks objects are constants that should be in capital letters
# pylint: disable=invalid-name

import sys
import os
import argparse

from wrapit.parser import add_doit_options
import pytest

from . import bias, cosegregation, matrix, call_windows, enrichment, radial_position, \
        permutation, plotting, pipeline, select_samples, slice, compaction


parser = argparse.ArgumentParser(
    prog='gamtools',
    description='GAMtools is a suite of tools for working with Genome'
    ' Architecture Mapping (GAM) data.')

subparsers = parser.add_subparsers(help='Please select a command:')
subparsers.required = True
subparsers.dest = 'command'

# Options for the 'bias' command

bias_parser = subparsers.add_parser(
    'bias',
    help='Calculate the bias in a set of matrices based on a genomic feature.')

bias_parser.add_argument(
    '-m', '--matrix-paths', metavar='MATRIX_PATH', nargs='+', required=True,
    help='Paths to matrix files for each chromosome')

bias_parser.add_argument(
    '-f', '--feature-path', metavar='FEATURE_PATH', required=True,
    help='Paths to bed file with one additional column giving the feature to \
    use for calculating biases')

bias_parser.add_argument(
    '-o', '--output-path', required=True,
    help='Path to save output matrix')

bias_parser.set_defaults(
    func=bias.calculate_bias_from_args)

# Options for the 'call_windows' command

call_windows_parser = subparsers.add_parser(
    'call_windows',
    help='Call positive windows for individual NPs')

call_windows_parser.add_argument(
    'coverage_file', type=argparse.FileType('r'),
    help='Input coverage file generated by bedtools multibamcov '
    '(or "-" to read from stdin)')

call_windows_parser.add_argument(
    '-d', '--details-file',
    help='If specified, write a table of fitting parameters to this path')

call_windows_parser.add_argument(
    '-o', '--output-file', type=argparse.FileType('w'),
    default=sys.stdout,
    help='Output segregation file to create '
    '(or "-" to write to stdout)')

# If a fixed threshold is specified, we don't need to plot the fits
seg_method = call_windows_parser.add_mutually_exclusive_group()
seg_method.add_argument(
    '-f', '--fitting-folder', metavar='FITTING_FOLDER',
    help='If specified, save the individual curve fittings to this folder')
seg_method.add_argument(
    '-x', '--fixed-threshold', metavar='NUM_READS',
    dest='fitting_function',
    type=call_windows.fixed_threshold_fitting_func,
    help='Instead of fitting each sample individually, use a fixed threshold '
    'to determine which windows are positive.')

call_windows_parser.set_defaults(
    func=call_windows.threshold_from_args,
    fitting_function=call_windows.signal_and_noise_fitting)

# Options for the 'compaction' command

compaction_parser = subparsers.add_parser(
    'compaction',
    help='Get the compaction of each genomic window from a segmentation file')


compaction_parser.add_argument(
    '-s', '--segregation-file', required=True,
    type=argparse.FileType('r'),
    help='A segregation file to use for calculating compaction'
    '(or - to read from stdin)')

compaction_parser.add_argument(
    '-o', '--output-file', required=True,
    type=argparse.FileType('w'),
    help='Output bedgraph file path (or - to write to stdout)')

compaction_parser.add_argument(
    '-n', '--no-blanks', default=False,
    action='store_true',
    help='Do not include lines with no value in the output.')

compaction_parser.set_defaults(func=compaction.compaction_from_args)



# Options for the 'convert' command

convert_parser = subparsers.add_parser(
    'convert',
    help='Convert between different GAM matrix formats')

convert_parser.add_argument(
    'input_file', metavar='INPUT_MATRIX',
    help='Input matrix file to convert')

convert_parser.add_argument(
    'output_file', metavar='OUTPUT_MATRIX',
    help='Output matrix file to write to '
    '(or "-" to write to stdout).')

convert_parser.add_argument(
    '-i', '--input-format',
    choices=matrix.INPUT_FORMATS,
    help='Input matrix file format (choose from: {})'.format(
        ', '.join(matrix.INPUT_FORMATS.keys())))

convert_parser.add_argument(
    '-o', '--output-format',
    choices=matrix.OUTPUT_FORMATS,
    help='Output matrix file format (choose from: {})'.format(
        ', '.join(matrix.OUTPUT_FORMATS.keys())))

convert_parser.add_argument(
    '-t', '--thresholds-file', metavar='THRESHOLDS_FILE',
    help='Thresholds file. If specified, any values lower than the specified'
    ' thresholds will be masked/excluded from the output file')

convert_parser.add_argument(
    '-w', '--windows-file', metavar='WINDOWS_FILE',
    help='File containing the genomic locations of matrix bins '
    '(only required if not specified in input matrix file).')

convert_parser.add_argument(
    '-r', '--region',
    help='Region covered by the input matrix '
    'required if -w/--windows-file is specified')

convert_parser.set_defaults(func=matrix.convert_from_args)


# Options for the 'enrichment' command

enrichment_parser = subparsers.add_parser(
    'enrichment',
    help='Calculate enrichments of SLICE interactions')

enrichment_parser.add_argument(
    '-i', '--interactions-file', required=True,
    help='A file containing all pairwise interactions between genomic windows')

enrichment_parser.add_argument(
    '-c', '--classes-file', required=True,
    help='A file containing classes assigned to each genomic window')

enrichment_parser.add_argument(
    '-o', '--output-prefix', default='enrichment_results',
    help='First part of the output file name')

# Either permuting the data or not permuting must be specified
perm_type = enrichment_parser.add_mutually_exclusive_group(required=True)
perm_type.add_argument(
    '-p', '--permutations', default=10, type=int, dest='num_permutations',
    help='Number of times to randomly permute the input file')
perm_type.add_argument(
    '-n',
    '--no-permute',
    action='store_const',
    const=0,
    dest='num_permutations',
    help='Do not permute the input file, instead calculate observed counts')

enrichment_parser.set_defaults(func=enrichment.enrichment_from_args)


# Options for the 'matrix' command

matrix_parser = subparsers.add_parser(
    'matrix',
    help='Generate a GAM matrix from a segregation file')

matrix_parser.add_argument(
    '-r', '--regions', metavar='REGION', required=True, nargs='+',
    help='Specific genomic regions to calculate matrices for. '
    'If one region is specified, a matrix is calculated for that region '
    'against itself. If more than one region is specified, a matrix is '
    'calculated for each region against the other. Regions are specified '
    'using UCSC browser syntax, i.e. "chr4" for the whole of chromosome '
    '4 or "chr4:100000-200000" for a sub-region of the chromosome.')

matrix_parser.add_argument(
    '-s', '--segregation_file', required=True,
    help='A segregation file to use as input')

matrix_parser.add_argument(
    '-f', '--output-format',
    choices=matrix.OUTPUT_FORMATS,
    help='Output matrix file format (choose from: {}, default is txt.gz)'.format(
        ', '.join(matrix.OUTPUT_FORMATS.keys())))

matrix_parser.add_argument(
    '-t',
    '--matrix-type',
    default='dprime',
    choices=cosegregation.MATRIX_TYPES,
    help='Method used to calculate the interaction matrix (choose from: '
    '{}, default is dprime)'.format(
        ', '.join(
            cosegregation.MATRIX_TYPES.keys())))

matrix_parser.add_argument(
    '-o', '--output-file', metavar='OUTPUT_FILE',
    help='Output matrix file. If not specified, new file will have the '
    'same name as the segregation file and an extension indicating the '
    'genomic region(s) and the matrix method')

matrix_parser.set_defaults(func=cosegregation.matrix_from_args)


# TODO: add Markus' matrix permutation code


# Options for 'permute_segregation' command

permute_seg_parser = subparsers.add_parser(
    'permute_segregation',
    help='Circularly permute the columns of a GAM segregation file')

permute_seg_parser.add_argument(
    '-s', '--segregation_file', required=True,
    type=argparse.FileType('r'),
    help='A segregation file to circularly permute '
    '(or - to read from stdin)')

permute_seg_parser.add_argument(
    '-o', '--output-file', required=True,
    type=argparse.FileType('w'),
    help='Output file path (or - to write to stdout)')

permute_seg_parser.set_defaults(func=permutation.permute_segregation_from_args)


# Options for 'plot_np' command

plot_np_parser = subparsers.add_parser(
    'plot_np',
    help='Plot the segregation results for a particular NP')

plot_np_parser.add_argument(
    '-w', '--bigwig-file', required=True,
    help='A bigwig file containing coverage information for the NP')
plot_np_parser.add_argument(
    '-b', '--bed-file', required=True,
    help='A bed file containing positive windows for the NP')
plot_np_parser.add_argument(
    '-g', '--genome-file', required=True,
    help='A file containing chromosome sizes for this genome')
plot_np_parser.add_argument(
    '-o', '--output-file', required=True, help='Output image file to create')

plot_np_parser.set_defaults(func=plotting.plot_np_from_args)


# Options for 'process' command

process_parser = subparsers.add_parser(
    'process_nps',
    help='Run a pipeline for mapping raw GAM sequencing data and '
    'calling positive windows')

process_parser.add_argument(
    'input_fastqs', metavar='INPUT_FASTQ', nargs='+',
    help='One or more input fastq files.')
process_parser.add_argument(
    '-g', '--genome-file', metavar='GENOME_FILE', required=True,
    help='File containing chromosome names and lengths')
process_parser.add_argument(
    '-o', '--output_dir', metavar='OUPUT_DIRECTORY',
    default=os.getcwd(),
    help='Write segregation, matrix etc. to this directory')
process_parser.add_argument(
    '-f', '--fittings_dir', metavar='FITTINGS_DIRECTORY',
    help='Write segregation curve fitting plots to this directory')
process_parser.add_argument(
    '-d', '--details-file',
    help='If specified, write a table of fitting parameters to this path')
process_parser.add_argument(
    '-i', '--bigwigs', action='append_const',
    dest='to_run', const='Converting bedgraph to bigwig',
    help='Make bigWig files.')
process_parser.add_argument(
    '-b', '--bigbeds', action='append_const',
    dest='to_run', const='Getting segregation bigBed',
    help='Make bed files of positive windows')
process_parser.add_argument(
    '-c', '--do-qc', action='append_const',
    dest='to_run', const='do_qc',
    help='Perform sample quality control.')
process_parser.add_argument(
    '-w', '--window-sizes', metavar='WINDOW_SIZE',
    default=[50000], type=int, nargs='+',
    help='One or more window sizes for calling positive windows')
process_parser.add_argument(
    '-m', '--matrix-sizes', metavar='MATRIX_SIZE',
    default=[], type=int, nargs='+',
    help='Resolutions for which linkage matrices should be produced.')
process_parser.add_argument(
    '--qc-window-size', type=int,
    help='Use this window size for qc (default is median window size).')
process_parser.add_argument(
    '--additional-qc-files', nargs='*', default=[],
    help='Any additional qc files to filter on')
process_parser.add_argument(
    '-q', '--minimum-mapq', metavar='MINIMUM_MAPQ', default=20, type=int,
    help='Filter out any mapped read with a mapping quality less than x '
    '(default is 20, use -q 0 for no filtering)')

# Add options for doit, the task runner engine (www.pydoit.org)
add_doit_options(process_parser,
                 ['dep_file', 'backend', 'verbosity',
                  'reporter', 'num_process', 'par_type',
                  'reset_dep'])


def get_script(script_name):
    """Get the absolute path to a script in the gamtools 'scripts' folder

    :param str script_name: Name of the script file.
    :returns: Absolute path to the script file
    """

    return os.path.join(os.path.dirname(__file__),
                        'data',
                        'scripts',
                        script_name)


def get_example(example_name):
    """Get the absolute path to a file in the gamtools 'examples' folder

    :param str example_name: Name of the file.
    :returns: Absolute path to the file
    """

    return os.path.join(os.path.dirname(__file__),
                        'data',
                        'examples',
                        example_name)

process_parser.set_defaults(func=pipeline.process_nps_from_args,
                            to_run=[
                                'Calling positive windows',
                                #'Filtering samples based on QC values',
                            ],
                            mapping_stats_script=get_script(
                                'mapping_stats.sh'),
                            example_parameters_file=get_example(
                                'qc_parameters.example.cfg'),
                            default_stats=['contamination_stats.txt', 'mapping_stats.txt',
                                           'quality_stats.txt', 'segregation_stats.txt'],
                            fitting_function=call_windows.signal_and_noise_fitting)

# Options for the 'radial_pos' command

radial_pos_parser = subparsers.add_parser(
    'radial_pos',
    help='Get the radial position of each genomic window from a segmentation file')


radial_pos_parser.add_argument(
    '-s', '--segregation-file', required=True,
    type=argparse.FileType('r'),
    help='A segregation file to use for calculating radial position'
    '(or - to read from stdin)')

radial_pos_parser.add_argument(
    '-o', '--output-file', required=True,
    type=argparse.FileType('w'),
    help='Output bedgraph file path (or - to write to stdout)')

radial_pos_parser.add_argument(
    '-n', '--no-blanks', default=False,
    action='store_true',
    help='Do not include lines with no value in the output.')


radial_pos_parser.set_defaults(func=radial_position.radial_position_from_args)

# Options for 'select' command

select_parser = subparsers.add_parser(
    'select',
    help='Select only certain samples from a segregation file')

select_parser.add_argument(
    '-s', '--segregation-file', required=True,
    help='A file containing the segregation of all samples')
select_parser.add_argument(
    '-d', '--drop-samples', action='store_true',
    help='Discard the listed samples (default: discard samples not in the list)')
select_parser.add_argument(
    '-n', '--sample-names', metavar='SAMPLE_NAME',
    required=True, nargs='*', help='Names of the samples to remove')
select_parser.add_argument(
    '-o', '--output-file', required=True,
    type=argparse.FileType('w'),
    help='Output file path (or - to write to stdout)')

select_parser.set_defaults(func=select_samples.select_samples_from_args)

# Options for 'slice' command

slice_parser = subparsers.add_parser(
    'slice',
    help='Use the SLICE library')

slice_parser.add_argument(
    '-s', '--segregation-file-path', required=True,
    help='Path to a segregation file')

slice_parser.add_argument(
    '-o', '--output-dir', required=True,
    help='Path to output directory')

slice_parser.add_argument(
    '-t', '--slice-thickness', default=0.22, type=float,
    help='Thickness of cryosections (same units as nuclear radius)')

slice_parser.add_argument(
    '-R', '--nuclear-radius', default=4.5, type=float,
    help='Radius of the nucleus (same units as slice thickness)')

slice_parser.set_defaults(func=slice.run_slice_from_args)
slice_parser.set_defaults(haploid_chroms=['chrX', 'chrY'])

# Options for 'test' command

test_parser = subparsers.add_parser(
    'test',
    help='Test that GAMtools is installed and configured correctly.')


test_directory = os.path.join(os.path.dirname(__file__), 'tests')

def test_function(args): #pylint: disable=unused-argument
    """Wrapper function to call py.test from arparse"""

    # TODO: Fix passing additional arguments to py.test
    test_args = sys.argv[2:]
    sys.exit(pytest.main([test_directory] + test_args))

test_parser.set_defaults(func=test_function)


def main():
    """Main gamtools function that is called by the gamtools command line script."""

    args = parser.parse_args()

    args.func(args)


if __name__ == '__main__':
    main()
