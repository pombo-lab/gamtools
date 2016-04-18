import os
from wrapit.api import run
import itertools
import inspect
from . import call_windows, segregation, qc, select_samples, cosegregation, utils

def swap_extension(input_fastq, extension):
    output_dir = os.path.dirname(input_fastq)
    basename, old_ext = os.path.splitext(os.path.basename(input_fastq))
    if old_ext == '.gz':
        basename, old_ext = os.path.splitext(basename)
    return os.path.join(output_dir, basename + extension)

def get_middle_value(_list):
    return sorted(_list)[len(_list)/2]

def pretty_resolution(window_size):
    return utils.format_genomic_distance(window_size, precision=0)

def coverage_path(base_folder, window_size):
    return os.path.join(base_folder, 'coverage_at_{0}.multibam'.format(pretty_resolution(window_size)))

def segregation_path(base_folder, window_size):
    return os.path.join(base_folder, 'segregation_at_{0}.multibam'.format(pretty_resolution(window_size)))

class input_file_mapping_tasks(object):

    def __init__(self, args):
        self.args = args

    def create_doit_tasks(self):

        tasks = []
        task_generators = []

        for possible_task in dir(self):
            if possible_task.startswith('task_'):

                task_obj = getattr(self, possible_task)()
                if inspect.isgenerator(task_obj):
                    task_generators.append(task_obj)
                else:
                    tasks.append(task_obj)

        task_generators.append((t for t in tasks))

        all_tasks = itertools.chain.from_iterable(task_generators)

        for task in all_tasks:
            yield task

    def task_map_fastq(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Mapping fastq',
                    "actions"  : ['bowtie2 -x genome -U %(dependencies)s | sed \'/XS : /d\' | samtools view -q %(minimum_mapq)s -F 4 -bS - > %(targets)s'],
                    "targets"  : [swap_extension(input_fastq, ".bam")],
                    "file_dep" : [input_fastq],
                  }

    def task_sort_bam(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Sorting bam',
                    "actions"  : ['samtools sort %(dependencies)s %(targets)s',
                                  'mv %(targets)s.bam %(targets)s'],
                    "targets"  : [swap_extension(input_fastq, ".sorted.bam")],
                    "file_dep" : [swap_extension(input_fastq, ".bam")],
                  }

    def task_remove_duplicates(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Removing duplications',
                    "actions"  : ['samtools rmdup -s %(dependencies)s %(targets)s',],
                    "targets"  : [swap_extension(input_fastq, ".rmdup.bam")],
                    "file_dep" : [swap_extension(input_fastq, ".sorted.bam")],
                  }

    def task_index_bam(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Indexing bam',
                    "actions"  : ['samtools index %(dependencies)s',],
                    "targets"  : [swap_extension(input_fastq, ".rmdup.bam.bai")],
                    "file_dep" : [swap_extension(input_fastq, ".rmdup.bam")],
                  }

    def task_make_bedgraph(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Converting bam to bedgraph',
                    "actions"  : ['genomeCoverageBed -bg -ibam %(dependencies)s -g %(genome_file)s > %(targets)s',
                                  'if [ -s "%(targets)s" ]; '
                                  'then create_empty_bedgraph %(genome_file)s %(targets)s; fi'],
                    "targets"  : [swap_extension(input_fastq, ".bedgraph")],
                    "file_dep" : [swap_extension(input_fastq, ".rmdup.bam")],
                  }

    def task_make_bigwig(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Converting bedgraph to bigwig',
                    "actions"  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                    "targets"  : [swap_extension(input_fastq, ".bw")],
                    "file_dep" : [swap_extension(input_fastq, ".bedgraph")],
                  }

    def task_make_output_dir(self):

        yield {
                "basename"     : 'Creating output directory',
                "actions"  : ['mkdir -p %(targets)s',],
                "targets"  : [self.args.output_dir],
                "uptodate" : [True],
              }

    def task_get_coverage(self):
        bamfiles = [ swap_extension(input_fastq, ".rmdup.bam") for input_fastq in self.args.input_fastqs ]
        for window_size in self.args.window_sizes:
            yield {
                    'basename' : 'Getting coverage',
                    'name'     :  '{0} windows'.format(pretty_resolution(window_size)),
                    'actions'  : ['echo "chrom start stop %(dependencies)s" > %(targets)s',
                                  'bash -i -c "multiBamCov -bams %(dependencies)s -bed <(bedtools makewindows -w {0} -g %(genome_file)s) >> %(targets)s"' .format(window_size)],
                    'targets'  : [coverage_path(self.args.output_dir, window_size)],
                    'file_dep' : bamfiles,
                    'task_dep' : [ 'Indexing bam', 'Creating output directory' ],
                   }

    def task_get_segregation(self):
        for window_size in self.args.window_sizes:

            input_coverage_file = coverage_path(self.args.output_dir, window_size)
            output_segregation_file = segregation_path(self.args.output_dir, window_size)

            task = {
                'basename' : 'Calling positive windows',
                'name'     : pretty_resolution(window_size),
                'targets'  : [output_segregation_file],
                'file_dep' : [input_coverage_file],
                'actions'  : [],
                }

            if self.args.fittings_dir:
                full_fitting_dir = os.path.join(self.args.output_dir, self.args.fittings_dir, '{0}'.format(pretty_resolution(window_size)))
                task['actions'].append('mkdir -pv {0}'.format(full_fitting_dir))
            else:
                full_fitting_dir = None

            task['actions'].append((call_windows.threshold_file, (
                                      input_coverage_file,
                                      output_segregation_file,
                                      full_fitting_dir,
                                      self.args.details_file,
                                      self.args.fitting_function)))

            yield task

    def task_get_segregation_bed(self):

        for window_size in self.args.window_sizes:

            input_segregation = segregation_path(self.args.output_dir, window_size)

            for input_fastq in self.args.input_fastqs:

                output_bed = swap_extension(input_fastq, ".segregation_{0}.bed".format(pretty_resolution(window_size)))

                yield {
                    'basename' : 'Getting segregation bed',
                    'name'     : '{0}, {1}'.format(pretty_resolution(window_size), input_fastq),
                    'actions'  : [(segregation.sample_segregation_to_bed,
                                   (input_segregation, swap_extension(input_fastq, '.rmdup.bam'), output_bed + '.unsorted')),
                                 'sort -k1,1 -k2,2n %(targets)s.unsorted > %(targets)s',
                                 'rm %(targets)s.unsorted',],
                    'targets'  : [output_bed],
                    'file_dep' : [segregation_path(self.args.output_dir, window_size)],
                    }

    def task_get_segregation_bigbed(self):

        for window_size in self.args.window_sizes:
            for input_fastq in self.args.input_fastqs:

                yield {
                    'basename' : 'Getting segregation bigBed',
                    'name'     : '{0}, {1}'.format(pretty_resolution(window_size), input_fastq),
                    # If the input bed is empty, bedToBigBed will fail. We can force it not to return an
                    # error code, but we must touch the target first to ensure it gets created.
                    'actions'  : ['touch %(targets)s',
                                  'bedToBigBed %(dependencies)s %(genome_file)s %(targets)s || true',],
                    'targets'  : [swap_extension(input_fastq, ".segregation_{0}.bb".format(pretty_resolution(window_size)))],
                    'file_dep' : [swap_extension(input_fastq, ".segregation_{0}.bed".format(pretty_resolution(window_size)))],
                      }

    def task_run_fastqc(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Running fastqc',
                    "actions"  : ['fastqc -q --extract %(dependencies)s'],
                    "targets"  : [qc.fastqc.fastqc_data_file(input_fastq)],
                    "file_dep" : [input_fastq],
                  }

    def task_run_fastq_screen(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Running fastq_screen',
                    "actions"  : ['fastq_screen --force --quiet %(dependencies)s'],
                    "targets"  : [qc.screen.screen_out_path(input_fastq)],
                    "file_dep" : [input_fastq],
                  }

    def task_extract_contamination(self):
        fastq_screen_files = [ qc.screen.screen_out_path(input_fastq) for input_fastq in self.args.input_fastqs ]
        return {
                'basename': 'Getting contamination stats',
                'actions' : [qc.screen.contamination_from_doit],
                'targets' : [os.path.join(self.args.output_dir, 'contamination_stats.txt')],
                'file_dep' : fastq_screen_files,
               }

    def task_extract_quality(self):
        fastqc_files = [qc.fastqc.fastqc_data_file(input_fastq) for input_fastq in self.args.input_fastqs ]
        return {
                'basename': 'Getting quality stats',
                'actions' : [qc.fastqc.quality_qc_from_doit],
                'targets' : [os.path.join(self.args.output_dir, 'quality_stats.txt')],
                'file_dep' : fastqc_files,
               }

    def task_mapping_stats(self):
        return {
                'basename': 'Getting mapping stats',
                'actions' : ['%(mapping_stats_script)s %(dependencies)s > %(targets)s'],
                'targets' : [os.path.join(self.args.output_dir, 'mapping_stats.txt')],
                'file_dep' : self.args.input_fastqs,
              }

    def task_segregation_stats(self):

        if self.args.qc_window_size is None:
            self.args.qc_window_size = get_middle_value(self.args.window_sizes)

        input_segregation = segregation_path(self.args.output_dir, self.args.qc_window_size)

        return {
                'basename': 'Getting segregation stats',
                'actions' : [qc.segregation.get_segregation_stats_doit],
                'targets' : [os.path.join(self.args.output_dir, 'segregation_stats.txt')],
                'file_dep': [input_segregation],
               }

    def task_merge_stats(self):

        files_to_merge = [os.path.join(self.args.output_dir, sf) for sf in self.args.default_stats] + self.args.additional_qc_files

        return {
                'basename': 'Merging stats files',
                'actions' : [qc.merge.merge_stats_from_doit],
                'targets' : [os.path.join(self.args.output_dir, 'merged_stats.txt')],
                'file_dep' : files_to_merge
               }

    def task_copy_qc_parameters(self):

        return {
                'basename': 'Creating QC parameters file with default values',
                'actions': ['cp %(example_parameters_file)s %(targets)s'],
                'targets': [os.path.join(self.args.output_dir, 'qc_parameters.cfg')],
                'uptodate' : [True],
               }

    def task_samples_passing_qc(self):

        return {
                'basename': 'Finding samples that pass QC',
                'actions': [qc.pass_qc.samples_passing_qc_from_doit],
                'targets': [os.path.join(self.args.output_dir, 'samples_passing_qc.txt')],
                'file_dep' : [os.path.join(self.args.output_dir, 'qc_parameters.cfg'),
                              os.path.join(self.args.output_dir, 'merged_stats.txt')],
               }

    def task_filter_segregation(self):

        passqc_file = os.path.join(self.args.output_dir, 'samples_passing_qc.txt')

        for window_size in self.args.window_sizes:

            segregation_file = segregation_path(self.args.output_dir, window_size)

            yield {
                'basename': 'Filtering samples based on QC values',
                'name'     : pretty_resolution(window_size),
                'targets'  : [swap_extension(segregation_file, '.passed_qc.multibam')],
                'file_dep' : [passqc_file, segregation_file],
                'actions'  : [select_samples.select_samples_from_doit]
                }

    def task_do_qc(self):
        return {'basename':'do_qc',
                'actions': None,
                'task_dep': ['Filtering samples based on QC values']}

    def task_linkage_matrices(self):

        #TODO: Allow specification of which chromosomes to generate matrices for
        chroms = ['chr{0}'.format(c) for c in range(1,20)]

        for window_size in self.args.matrix_sizes:

            for chrom in chroms:

                segregation_file = segregation_path(self.args.output_dir, window_size)
                if 'do_qc' in self.args.to_run:
                    segregation_file = swap_extension(segregation_file, '.passed_qc.multibam')

                matrix_base = '{seg_file}.dprime.{chrom}.txt.gz'.format(
                                                  seg_file=os.path.basename(segregation_file),
                                                  chrom=chrom)
                matrix_out = os.path.join(self.args.output_dir, 'dprime_matrices', '{}'.format(pretty_resolution(window_size)), matrix_base)

                yield {'basename' : 'Calculating linkage matrix',
                       'name'     : '{0} at {1}'.format(chrom, pretty_resolution(window_size)),
                       'targets'  : [matrix_out],
                       'file_dep' : [segregation_file],
                       'actions' : ['mkdir -p $(dirname %(targets)s)',
                                    (cosegregation.matrix_from_doit, (matrix_out, segregation_file, [chrom]))],}


def process_nps_from_args(args):

    task_dict = {
        'input_file_tasks': input_file_mapping_tasks(args)
    }
    run(task_dict, args, args.to_run)

