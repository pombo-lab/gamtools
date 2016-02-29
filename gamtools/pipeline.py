import os
from wrapit.api import run
import itertools
import inspect
from . import call_windows, segmentation, qc

def with_extension(input_fastq, extension):
    output_dir = os.path.dirname(input_fastq)
    basename, old_ext = os.path.splitext(os.path.basename(input_fastq))
    if old_ext == '.gz':
        basename, old_ext = os.path.splitext(basename)
    return os.path.join(output_dir, basename + extension)

def segmentation_path(base_folder, window_size):

    return os.path.join(base_folder, 'segmentation_at_{0}bp.multibam'.format(window_size))

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
                    "targets"  : [with_extension(input_fastq, ".bam")],
                    "file_dep" : [input_fastq],
                  }

    def task_sort_bam(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Sorting bam',
                    "actions"  : ['samtools sort %(dependencies)s %(targets)s',
                                  'mv %(targets)s.bam %(targets)s'],
                    "targets"  : [with_extension(input_fastq, ".sorted.bam")],
                    "file_dep" : [with_extension(input_fastq, ".bam")],
                  }

    def task_remove_duplicates(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Removing duplications',
                    "actions"  : ['samtools rmdup -s %(dependencies)s %(targets)s',],
                    "targets"  : [with_extension(input_fastq, ".rmdup.bam")],
                    "file_dep" : [with_extension(input_fastq, ".sorted.bam")],
                  }

    def task_index_bam(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Indexing bam',
                    "actions"  : ['samtools index %(dependencies)s',],
                    "targets"  : [with_extension(input_fastq, ".rmdup.bam.bai")],
                    "file_dep" : [with_extension(input_fastq, ".rmdup.bam")],
                  }

    def task_make_bedgraph(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Converting bam to bedgraph',
                    "actions"  : ['genomeCoverageBed -bg -ibam %(dependencies)s -g %(genome_file)s > %(targets)s',
                                  'if [ -s "%(targets)s" ]; '
                                  'then create_empty_bedgraph %(genome_file)s %(targets)s; fi'],
                    "targets"  : [with_extension(input_fastq, ".bedgraph")],
                    "file_dep" : [with_extension(input_fastq, ".rmdup.bam")],
                  }

    def task_make_bigwig(self):
        for input_fastq in self.args.input_fastqs:
            yield {
                    "name"     : input_fastq,
                    "basename" : 'Converting bedgraph to bigwig',
                    "actions"  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                    "targets"  : [with_extension(input_fastq, ".bw")],
                    "file_dep" : [with_extension(input_fastq, ".bedgraph")],
                  }

    def task_make_output_dir(self):

        yield {
                "basename"     : 'Creating output directory',
                "actions"  : ['mkdir -pv %(output_dir)s',],
                "targets"  : ['%(output_dir)s'],
              }

    def task_get_coverage(self):
        bamfiles = [ with_extension(input_fastq, ".rmdup.bam") for input_fastq in self.args.input_fastqs ]
        for window_size in self.args.window_sizes:
            yield {
                    'basename' : 'Getting coverage',
                    'name'     :  '{0}bp windows'.format(window_size),
                    'actions'  : ['echo "chrom start stop %(dependencies)s" > %(targets)s',
                                  'bash -i -c "multiBamCov -bams %(dependencies)s -bed <(bedtools makewindows -w {0} -g %(genome_file)s) >> %(targets)s"' .format(window_size)],
                    'targets'  : [os.path.join(self.args.output_dir, 'coverage_at_{0}bp.multibam'.format(window_size))],
                    'file_dep' : bamfiles,
                    'task_dep' : [ 'Indexing bam', 'Creating output directory' ],
                   }

    def task_get_segmentation(self):
        for window_size in self.args.window_sizes:

            input_coverage_file = os.path.join(self.args.output_dir, 'coverage_at_{0}bp.multibam'.format(window_size))
            output_segmentation_file = segmentation_path(self.args.output_dir, window_size)

            task = {
                'basename' : 'Calling positive windows',
                'name'     : '{0}bp'.format(window_size),
                'targets'  : [output_segmentation_file],
                'file_dep' : [input_coverage_file],
                'actions'  : [],
                }

            if self.args.fittings_dir:
                full_fitting_dir = os.path.join(self.args.output_dir, self.args.fittings_dir, '{0}bp'.format(window_size))
                task['actions'].append('mkdir -pv {0}'.format(full_fitting_dir))
            else:
                full_fitting_dir = None

            task['actions'].append((call_windows.threshold_file, (
                                      input_coverage_file,
                                      output_segmentation_file,
                                      full_fitting_dir,
                                      self.args.details_file,
                                      self.args.fitting_function)))

            yield task

    def task_get_segmentation_bed(self):

        for window_size in self.args.window_sizes:

            input_segmentation = segmentation_path(self.args.output_dir, window_size)

            for input_fastq in self.args.input_fastqs:

                output_bed = with_extension(input_fastq, ".segmentation_{0}bp.bed".format(window_size))

                yield {
                    'basename' : 'Getting segmentation bed',
                    'name'     : '{0}bp, {1}'.format(window_size, input_fastq),
                    'actions'  : [(segmentation.sample_segmentation_to_bed,
                                   (input_segmentation, with_extension(input_fastq, '.rmdup.bam'), output_bed + '.unsorted')),
                                 'sort -k1,1 -k2,2n %(targets)s.unsorted > %(targets)s',
                                 'rm %(targets)s.unsorted',],
                    'targets'  : [output_bed],
                    'file_dep' : [os.path.join(self.args.output_dir, 'segmentation_at_{0}bp.multibam'.format(window_size))],
                    }

    def task_get_segmentation_bigbed(self):

        for window_size in self.args.window_sizes:
            for input_fastq in self.args.input_fastqs:

                yield {
                    'basename' : 'Getting segmentation bigBed',
                    'name'     : '{0}bp, {1}'.format(window_size, input_fastq),
                    # If the input bed is empty, bedToBigBed will fail. We can force it not to return an
                    # error code, but we must touch the target first to ensure it gets created.
                    'actions'  : ['touch %(targets)s',
                                  'bedToBigBed %(dependencies)s %(genome_file)s %(targets)s || true',],
                    'targets'  : [with_extension(input_fastq, ".segmentation_{0}bp.bb".format(window_size))],
                    'file_dep' : [with_extension(input_fastq, ".segmentation_{0}bp.bed".format(window_size))],
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

    def task_segmentation_stats(self):

        input_segmentation = segmentation_path(self.args.output_dir, self.args.qc_window_size)

        return {
                'basename': 'Getting segmentation stats',
                'actions' : [qc.segmentation.get_segmentation_stats_doit],
                'targets' : [os.path.join(self.args.output_dir, 'segmentation_stats.txt')],
                'file_dep': [input_segmentation]
               }

    def task_merge_stats(args):

        files_to_merge = [os.path.join(args.output_dir, sf) for sf in args.default_stats] + args.additional_qc_files

        return {
                'basename': 'Merging stats files',
                'actions' : [qc.merge.merge_stats_from_doit],
                'targets' : [os.path.join(args.output_dir, 'merged_stats.txt')],
                'file_dep' : files_to_merge
               }

def process_nps_from_args(args):

    #from doit.loader import load_tasks
    #print load_tasks()
    task_dict = {
        'input_file_tasks': input_file_mapping_tasks(args)
    }
    #run(task_dict, args, ['list'])
    run(task_dict, args, args.to_run)

