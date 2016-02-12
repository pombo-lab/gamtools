

import argparse

parser = argparse.ArgumentParser(
    description='GAMtools is a suite of tools for working with Genome'
    ' Architecture Mapping (GAM) data.')

subparsers = parser.add_subparsers(help='Please select a command:')

call_windows_parser = subparsers.add_parser(
    'call_windows',
    help='Call positive windows for individual NPs')

convert_parser = subparsers.add_parser(
    'convert',
    help='Convert between different GAM matrix formats')

enrichment_parser = subparsers.add_parser(
    'enrichment',
    help='Calculate enrichments of SLICE interactions')

matrix_parser = subparsers.add_parser(
    'matrix',
    help='Generate a GAM matrix from a segmentation file')

permute_matrix_parser = subparsers.add_parser(
    'permute_matrix',
    help='Circularly permute a GAM matrix')

permute_seg_parser = subparsers.add_parser(
    'permute_segmentation',
    help='Circularly permute a GAM segmentation file')

plot_np_parser = subparsers.add_parser(
    'plot_np',
    help='Plot the segmentation results for a particular NP')

process_parser = subparsers.add_parser(
    'process',
    help='Run a pipeline for processing raw GAM sequencing data')

seg_stats_parser = subparsers.add_parser(
    'segmentation_stats',
    help='Calculate QC statistics from a segmentation file')

select_parser = subparsers.add_parser(
    'select',
    help='Select only certain samples from a segmentation file')

threshold_parser = subparsers.add_parser(
    'threshold_matrix',
    help='Apply a threshold to a SLICE interaction matrix')

if __name__ == '__main__':
    parser.parse_args()
