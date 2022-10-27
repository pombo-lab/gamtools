"""
===================
The pipeline module
===================

The pipeline module handles the :docs:`tools/process_nps` command.

The process_nps command uses doit_ as a task engine, to determine
what actions (or tasks) need to be done and in what order. Since
the pipeline is part of a larger command line program, it needs some
way to take arguments/parameters outside of the normal doit framework.
This functionality is provided by wrapit_ - which is an extension of
doit that allows additional parameters to be passed to tasks.

:class:`InputFileMappingTasks` is a class that saves the input
parameters and uses them to create tasks. Any method of the class
that begins with `task_` is considered a task and will be passed
to doit.

.. _doit: http://pydoit.org/

"""

import os
import sys
import subprocess
import re
import itertools
import inspect

from distutils.version import LooseVersion #pylint: disable=no-name-in-module,import-error
from wrapit.api import run

from . import call_windows, segregation, qc, select_samples, cosegregation, utils


def swap_extension(input_path, extension):
    """Replace the file extension of a file path.

    :param str input_path: Input file path
    :param str extension: New extension to use
    :returns: Input file path with the extension replaced.
    :rtype: str
    """

    output_dir = os.path.dirname(input_path)
    basename, old_ext = os.path.splitext(os.path.basename(input_path))
    if old_ext == '.gz':
        basename, old_ext = os.path.splitext(basename)
    return os.path.join(output_dir, basename + extension)


def get_middle_value(my_list):
    """Return the middle value from a list after sorting.

    :param list my_list: List of sortable values"""

    return sorted(my_list)[len(my_list) // 2]


def pretty_resolution(window_size):
    """Convert an integer resolution in base-pairs to a nicely formatted string.

    :param int window_size: Integer resolution in base pairs.
    :returns: Formatted resolution."""

    return utils.format_genomic_distance(window_size, precision=0)


def bp_coverage_path(base_folder, window_size):
    """Get the path to a bp coverage table given a base folder and resolution.

    :param str base_folder: Path to the folder containing the bp coverage table.
    :param int window_size: Resolution in base pairs
    :returns: Path to bp coverage table
    """

    return os.path.join(
        base_folder,
        'bp_coverage_at_{0}.table'.format(
            pretty_resolution(window_size)))


def individual_bp_coverage_path(window_size, bam_path):
    """Get the path to a bp coverage table for one specific bam

    :param int window_size: Resolution in base pairs
    :param str bam_path: Path to input bam for bp coverage calculation
    :returns: Path to bp coverage table
    """

    new_suffix = '.bp_at_{}.bed'.format(
        pretty_resolution(window_size))
    return swap_extension(bam_path, new_suffix)


def segregation_path(base_folder, window_size):
    """Get the path to a segregation table given a base folder and resolution.

    :param str base_folder: Path to the folder containing the segregation table.
    :param int window_size: Resolution in base pairs
    :returns: Path to segregation table
    """

    return os.path.join(
        base_folder,
        'segregation_at_{0}.table'.format(
            pretty_resolution(window_size)))


def get_samtools_version():
    """Get the version number for the installed version of samtools"""

    try:

        with subprocess.Popen('samtools',
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE) as proc:

            output, err = proc.communicate() #pylint: disable=unused-variable
            proc.wait()

    except OSError:
        sys.exit('samtools is either not installed or not present in $PATH. '
                 'Please install samtools to continue.')

    vnum_regexp = b'^Version:\s(?P<version>(\d+\.){1,2}(\d+))'

    regexp_matches = re.search(vnum_regexp, err, re.MULTILINE)

    return regexp_matches.group('version').decode('UTF-8')


def get_samtools_sort_actions():
    """Return the appropriate doit function for using the installed version of samtools sort"""

    samtools_vnum = get_samtools_version()
    samtools_vnum = LooseVersion(samtools_vnum)

    if samtools_vnum < LooseVersion('0.1.19'):
        samtools_actions = ['samtools sort %(dependencies)s %(targets)s.tmp',
                            'mv %(targets)s.tmp.bam %(targets)s']

    elif samtools_vnum < LooseVersion('1.0'):
        samtools_actions = ['samtools sort -f %(dependencies)s %(targets)s']

    else:
        samtools_actions = [
            'samtools sort -o %(targets)s -T %(targets)s.tmp %(dependencies)s']

    return samtools_actions


class InputFileMappingTasks():
    """Class for generating doit tasks from command-line arguments.

    GAMtools "process_nps" command generates a set of doit tasks at
    runtime based on a set of parameters passed via the command-line.
    These parameters are handled by argparse (see main.py for the
    parser definitions), yielding an object called args. Tasks
    need access to args at runtime, so they must be generated
    via a class, such that args is always accessible via self.args.

    :param args: An argparse Namespace object containing parameters \
            passed via the command line
    """

    def __init__(self, args):
        self.args = args

    def create_doit_tasks(self):
        """Generator function that yields doit tasks."""

        tasks = []
        task_generators = []

        for possible_task in dir(self):

            # Only attributes starting with task_ are task generators
            if possible_task.startswith('task_'):

                # Call each function starting with task_
                task_obj = getattr(self, possible_task)()

                # We cannot mix generators and task dictionaries
                # in the same iterable
                if inspect.isgenerator(task_obj):
                    task_generators.append(task_obj)
                else:
                    tasks.append(task_obj)

        # Aggregate task dictionaries into a generator and add them
        # to task_generators
        task_generators.append((t for t in tasks))

        # Turn the list of generators into a single iterable
        all_tasks = itertools.chain.from_iterable(task_generators)

        for task in all_tasks:
            yield task

    def task_map_fastq(self):
        """Task to map fastq files using bowtie."""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Mapping fastq',
                "actions": ['bowtie2 -x genome -U %(dependencies)s | sed \'/XS : /d\' | '
                            'samtools view -q %(minimum_mapq)s -F 4 -bS - > %(targets)s'],
                "targets": [swap_extension(input_fastq, ".bam")],
                "file_dep": [input_fastq],
            }

    def task_sort_bam(self):
        """Task to sort bam files"""

        sort_actions = get_samtools_sort_actions()

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Sorting bam',
                "actions": sort_actions,
                "targets": [swap_extension(input_fastq, ".sorted.bam")],
                "file_dep": [swap_extension(input_fastq, ".bam")],
            }

    def task_remove_duplicates(self):
        """Task to remove PCR duplicates from sorted bam files"""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Removing duplications',
                "actions": ['samtools rmdup -s %(dependencies)s %(targets)s', ],
                "targets": [swap_extension(input_fastq, ".rmdup.bam")],
                "file_dep": [swap_extension(input_fastq, ".sorted.bam")],
            }

    def task_index_bam(self):
        """Task to generate indexes for sorted bam files"""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Indexing bam',
                "actions": ['samtools index %(dependencies)s', ],
                "targets": [swap_extension(input_fastq, ".rmdup.bam.bai")],
                "file_dep": [swap_extension(input_fastq, ".rmdup.bam")],
            }

    def task_make_bedgraph(self):
        """Task to create bedgraphs from indexed and sorted bam files"""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Converting bam to bedgraph',
                "actions": ['genomeCoverageBed -bg -ibam %(dependencies)s '
                            '-g %(genome_file)s > %(targets)s',
                            'if [ ! -s "%(targets)s" ]; '
                            'then echo "Creating empty bedgraph %(targets)s";'
                            'create_empty_bedgraph %(genome_file)s %(targets)s; fi'],
                "targets": [swap_extension(input_fastq, ".bedgraph")],
                "file_dep": [swap_extension(input_fastq, ".rmdup.bam")],
            }

    def task_make_bigwig(self):
        """Task to create bigwigs from bedgraphs"""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Converting bedgraph to bigwig',
                "actions": ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                "targets": [swap_extension(input_fastq, ".bw")],
                "file_dep": [swap_extension(input_fastq, ".bedgraph")],
            }

    def task_make_output_dir(self):
        """Task to create a directory"""

        yield {
            "basename": 'Creating output directory',
            "actions": ['mkdir -p %(targets)s', ],
            "targets": [self.args.output_dir],
            "uptodate": [True],
        }

    def task_get_bp_coverage(self):
        """Task to get read coverage per window for a range of resolutions"""

        for input_fastq in self.args.input_fastqs:
            for window_size in self.args.window_sizes:
                input_bam = swap_extension(input_fastq, ".rmdup.bam")
                output_cov = individual_bp_coverage_path(window_size, input_bam)
                yield {
                    'basename': 'Getting bp coverage',
                    'name': output_cov,
                    'actions': ['echo "chrom start stop %(dependencies)s" > %(targets)s',
                                'bash -i -c "bedtools coverage -b %(dependencies)s '
                                '-a <(bedtools makewindows -w {0} -g %(genome_file)s) '
                                '| cut -d \'\t\' -f 1,2,3,7 '
                                '>> %(targets)s"' .format(window_size)],
                    'targets': [output_cov],
                    'file_dep': [input_bam],
                    'task_dep': ['Indexing bam', 'Creating output directory'],
                }


    def task_merge_bp_coverage(self):
        """Task to merge bp coverage files for individual bams into a bp coverage table"""

        input_bamfiles = [swap_extension(input_fastq, ".rmdup.bam")
                          for input_fastq in self.args.input_fastqs]

        for window_size in self.args.window_sizes:

            input_cov_files = [individual_bp_coverage_path(window_size, input_bam)
                               for input_bam in input_bamfiles]

            yield {
                'basename': 'Merging bp coverage',
                'name': '{0} windows'.format(pretty_resolution(window_size)),
                'actions': [call_windows.merge_reads_coverage],
                'targets': [bp_coverage_path(self.args.output_dir, window_size)],
                'file_dep': input_cov_files,
            }


    def task_get_segregation(self):
        """Task to call positive windows from coverage tables for a range of resolutions"""

        for window_size in self.args.window_sizes:

            input_coverage_file = bp_coverage_path(
                self.args.output_dir, window_size)
            output_segregation_file = segregation_path(
                self.args.output_dir, window_size)

            task = {
                'basename': 'Calling positive windows',
                'name': pretty_resolution(window_size),
                'targets': [output_segregation_file],
                'file_dep': [input_coverage_file],
                'actions': [],
            }

            if self.args.fitting_folder:
                full_fitting_dir = os.path.join(
                    self.args.output_dir, self.args.fitting_folder, '{0}'.format(
                        pretty_resolution(window_size)))
                task['actions'].append(
                    'mkdir -pv {0}'.format(full_fitting_dir))
            else:
                full_fitting_dir = None

            task['actions'].append((call_windows.threshold_file, (
                input_coverage_file,
                output_segregation_file,
                full_fitting_dir,
                self.args.details_file,
                self.args.fitting_function,
                self.args.min_read_threshold)))

            yield task

    def task_get_segregation_bed(self):
        """Task to get bed files of positive windows for each NP at each resolution"""

        for window_size in self.args.window_sizes:

            input_segregation = segregation_path(
                self.args.output_dir, window_size)

            for input_fastq in self.args.input_fastqs:

                output_bed = swap_extension(
                    input_fastq, ".segregation_{0}.bed".format(
                        pretty_resolution(window_size)))

                yield {
                    'basename': 'Getting segregation bed',
                    'name': '{0}, {1}'.format(pretty_resolution(window_size), input_fastq),
                    'actions': [(segregation.sample_segregation_to_bed,
                                 (input_segregation,
                                  swap_extension(input_fastq, '.rmdup.bam'),
                                  output_bed + '.unsorted')),
                                'sort -k1,1 -k2,2n %(targets)s.unsorted > %(targets)s',
                                'rm %(targets)s.unsorted', ],
                    'targets': [output_bed],
                    'file_dep': [segregation_path(self.args.output_dir, window_size)],
                }

    def task_get_segregation_bigbed(self):
        """Task to get bigbed files of positive windows for each NP at each resolution"""

        for window_size in self.args.window_sizes:
            for input_fastq in self.args.input_fastqs:

                yield {
                    'basename': 'Getting segregation bigBed',
                    'name': '{0}, {1}'.format(pretty_resolution(window_size), input_fastq),
                    # If the input bed is empty, bedToBigBed will fail. We can force it not
                    # to return an error code, but we must touch the target first to ensure
                    # it gets created.
                    'actions': ['touch %(targets)s',
                                'bedToBigBed %(dependencies)s %(genome_file)s %(targets)s || '
                                'true', ],
                    'targets': [swap_extension(input_fastq,
                                               ".segregation_{0}.bb".format(
                                                   pretty_resolution(window_size)))],
                    'file_dep': [swap_extension(input_fastq,
                                                ".segregation_{0}.bed".format(
                                                    pretty_resolution(window_size)))],
                }

    def task_run_fastqc(self):
        """Task to run fastqc for each input fastq"""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Running fastqc',
                "actions": ['fastqc -q --extract %(dependencies)s'],
                "targets": [qc.fastqc.fastqc_data_file(input_fastq)],
                "file_dep": [input_fastq],
            }

    def task_run_fastq_screen(self):
        """Task to run fastq_screen for each input fastq"""

        for input_fastq in self.args.input_fastqs:
            yield {
                "name": input_fastq,
                "basename": 'Running fastq_screen',
                "actions": ['fastq_screen --force --quiet %(dependencies)s'],
                "targets": [qc.screen.screen_out_path(input_fastq)],
                "file_dep": [input_fastq],
            }

    def task_extract_contamination(self):
        """Task to summarize results of fastq_screen runs"""

        fastq_screen_files = [qc.screen.screen_out_path(
            input_fastq) for input_fastq in self.args.input_fastqs]
        return {
            'basename': 'Getting contamination stats',
            'actions': [
                qc.screen.contamination_from_doit],
            'targets': [
                os.path.join(
                    self.args.output_dir,
                    'contamination_stats.txt')],
            'file_dep': fastq_screen_files,
        }

    def task_extract_quality(self):
        """Task to summarize results of fastqc runs"""

        fastqc_files = [qc.fastqc.fastqc_data_file(
            input_fastq) for input_fastq in self.args.input_fastqs]
        return {
            'basename': 'Getting quality stats',
            'actions': [
                qc.fastqc.quality_qc_from_doit],
            'targets': [
                os.path.join(
                    self.args.output_dir,
                    'quality_stats.txt')],
            'file_dep': fastqc_files,
        }

    def task_mapping_stats(self):
        """Task to count numbers of sequenced/mapped/duplicated reads"""

        return {
            'basename': 'Getting mapping stats',
            'actions': ['%(mapping_stats_script)s %(dependencies)s > %(targets)s'],
            'targets': [
                os.path.join(
                    self.args.output_dir,
                    'mapping_stats.txt')],
            'file_dep': self.args.input_fastqs,
        }

    def task_segregation_stats(self):
        """Task to calculate summary statistics for positive window calling"""

        if self.args.qc_window_size is None:
            self.args.qc_window_size = get_middle_value(self.args.window_sizes)

        input_segregation = segregation_path(
            self.args.output_dir, self.args.qc_window_size)

        return {
            'basename': 'Getting segregation stats',
            'actions': [
                qc.segregation.get_segregation_stats_doit],
            'targets': [
                os.path.join(
                    self.args.output_dir,
                    'segregation_stats.txt')],
            'file_dep': [input_segregation],
        }

    def task_merge_stats(self):
        """Task to make a single QC statistics table with one row per NP"""

        files_to_merge = [
            os.path.join(
                self.args.output_dir,
                sf) for sf in self.args.default_stats] + self.args.additional_qc_files

        return {
            'basename': 'Merging stats files',
            'actions': [
                qc.merge.merge_stats_from_doit],
            'targets': [
                os.path.join(
                    self.args.output_dir,
                    'merged_stats.txt')],
            'file_dep': files_to_merge}

    def task_copy_qc_parameters(self):
        """Task to make a qc_parameters.cfg file if it doesn't exist"""

        return {
            'basename': 'Creating QC parameters file with default values',
            'actions': ['cp %(example_parameters_file)s %(targets)s'],
            'targets': [
                os.path.join(
                    self.args.output_dir,
                    'qc_parameters.cfg')],
            'uptodate': [True],
        }

    def task_samples_passing_qc(self):
        """Task to identify NPs passing/failing QC"""

        return {
            'basename': 'Finding samples that pass QC',
            'actions': [qc.pass_qc.samples_passing_qc_from_doit],
            'targets': [os.path.join(self.args.output_dir, 'samples_passing_qc.txt')],
            'file_dep': [
                os.path.join(self.args.output_dir, 'qc_parameters.cfg'),
                os.path.join(self.args.output_dir, 'merged_stats.txt')], }

    def task_filter_segregation(self):
        """Task to exclude NPs failing QC from segregation tables"""

        passqc_file = os.path.join(
            self.args.output_dir,
            'samples_passing_qc.txt')

        for window_size in self.args.window_sizes:

            segregation_file = segregation_path(
                self.args.output_dir, window_size)

            yield {
                'basename': 'Filtering samples based on QC values',
                'name': pretty_resolution(window_size),
                'targets': [swap_extension(segregation_file, '.passed_qc.table')],
                'file_dep': [passqc_file, segregation_file],
                'actions': [select_samples.select_samples_from_doit]
            }

    def task_do_qc(self):
        """Empty target to force all QC tasks to run"""

        return {'basename': 'do_qc',
                'actions': None,
                'task_dep': ['Filtering samples based on QC values']}

    def task_linkage_matrices(self):
        """Task to calculate linkage matrices at specified resolutions"""

        # TODO: Allow specification of which chromosomes to generate matrices
        # for
        chroms = ['chr{0}'.format(c) for c in range(1, 20)]

        for window_size in self.args.matrix_sizes:

            for chrom in chroms:

                segregation_file = segregation_path(
                    self.args.output_dir, window_size)
                if 'do_qc' in self.args.to_run:
                    segregation_file = swap_extension(
                        segregation_file, '.passed_qc.table')

                matrix_base = '{seg_file}.dprime.{chrom}.txt.gz'.format(
                    seg_file=os.path.basename(segregation_file),
                    chrom=chrom)
                matrix_out = os.path.join(
                    self.args.output_dir, 'dprime_matrices', '{}'.format(
                        pretty_resolution(window_size)), matrix_base)

                yield {'basename': 'Calculating linkage matrix',
                       'name': '{0} at {1}'.format(chrom, pretty_resolution(window_size)),
                       'targets': [matrix_out],
                       'file_dep': [segregation_file],
                       'actions': ['mkdir -p $(dirname %(targets)s)',
                                   (cosegregation.matrix_from_doit,
                                    (matrix_out, segregation_file, [chrom]))], }

def check_resolution_consistency(args):
    """
    Ensure that matrix_sizes and qc_window_size are subsets of window_sizes.
    """

    all_window_sizes = list(args.window_sizes)
    all_window_sizes.extend(args.matrix_sizes)
    if args.qc_window_size is not None:
        all_window_sizes.append(args.qc_window_size)

    args.window_sizes = sorted(list(set(all_window_sizes)))

def process_nps_from_args(args):
    """Wrapper function to call doit from argparse"""

    check_resolution_consistency(args)

    if args.matrix_sizes:
        args.to_run.append('Calculating linkage matrix')

    task_dict = {
        'input_file_tasks': InputFileMappingTasks(args)
    }
    run(task_dict, args, args.to_run)
