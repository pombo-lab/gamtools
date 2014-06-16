#! /usr/bin/env python

from doit.loader import generate_tasks
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
import argparse
import os
import sys


parser = argparse.ArgumentParser(description='Map GAM-seq reads and create BigWigs and fastQC files')
parser.add_argument('input_fastqs', metavar='INPUT_FASTQ', nargs='+', help='One or more input fastq files.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', help='File containing chromosome names and lengths')
parser.add_argument('-p','--num_process', metavar='NUM_PROCESSES', type=int, default=1, help='Number of tasks to execute in parallel')
parser.add_argument('-o','--output_dir', metavar='OUPUT_DIRECTORY', help='Write segmentation, matrix etc. to this directory')
parser.add_argument('-f','--fittings_dir', metavar='FITTINGS_DIR', help='Write segmentation curve fitting plots to this directory')
parser.add_argument('-b','--bigbeds', action='store_true', help='Make bed files for segmentation')
parser.add_argument('-w','--window_sizes', metavar='WINDOW_SIZE', default=[1000,5000,10000,50000,100000,500000], type=int, nargs='+', help='File containing chromosome names and lengths')
parser.add_argument('-q','--minimum-mapq', metavar='MINIMUM_MAPQ', default=0, type=int, help='Filter out any mapped read with a mapping quality less than x (default is not to filter based on quality')

args = parser.parse_args()

def with_extension(input_fastq, extension): 
    output_dir = os.path.dirname(input_fastq)
    basename = os.path.splitext(os.path.basename(input_fastq))[0]
    return os.path.join(output_dir, basename + extension)

def map_fastq(input_fastqs):
    for input_fastq in input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Mapping fastq',
                "actions"  : ['bowtie2 -x genome -U %(dependencies)s | sed \'/XS : /d\' | samtools view -q %(minimum_mapq)s -F 4 -bS - > %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".bam")],
                "file_dep" : [input_fastq],
              }

def sort_bam(input_fastqs):
    for input_fastq in input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Sorting bam',
                "actions"  : ['samtools sort %(dependencies)s %(targets)s',
                              'mv %(targets)s.bam %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".sorted.bam")],
                "file_dep" : [with_extension(input_fastq, ".bam")],
              }

def remove_duplicates(input_fastqs):
    for input_fastq in input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Removing duplications',
                "actions"  : ['samtools rmdup -s %(dependencies)s %(targets)s',],
                "targets"  : [with_extension(input_fastq, ".rmdup.bam")],
                "file_dep" : [with_extension(input_fastq, ".sorted.bam")],
              }

def index_bam(input_fastqs):
    for input_fastq in input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Indexing bam',
                "actions"  : ['samtools index %(dependencies)s',],
                "targets"  : [with_extension(input_fastq, ".rmdup.bam.bai")],
                "file_dep" : [with_extension(input_fastq, ".rmdup.bam")],
              }

def make_bedgraph(input_fastqs):
    for input_fastq in input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Converting bam to bedgraph',
                "actions"  : ['genomeCoverageBed -bg -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".bedgraph")],
                "file_dep" : [with_extension(input_fastq, ".rmdup.bam")],
              }

def make_bigwig(input_fastqs):
    for input_fastq in input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Converting bedgraph to bigwig',
                "actions"  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".bw")],
                "file_dep" : [with_extension(input_fastq, ".bedgraph")],
              }

def get_coverage(args):
    bamfiles = [ with_extension(input_fastq, ".rmdup.bam") for input_fastq in args.input_fastqs ]
    for window_size in args.window_sizes:
        yield {
                'basename' : 'Getting coverage',
                'name' :  '{0}bp windows'.format(window_size),
                'actions' : ['echo "chrom start stop %(dependencies)s" > %(targets)s',
                             'bash -i -c "multiBamCov -bams %(dependencies)s -bed <(bedtools makewindows -w {0} -g %(genome_file)s) >> %(targets)s"' .format(window_size)],
                'targets' : [os.path.join(args.output_dir, 'coverage_at_{0}bp.multibam'.format(window_size))],
                'file_dep' : bamfiles,
                'task_dep' : [ 'Indexing bam' ],
               }

def get_segmentation(args):
    for window_size in args.window_sizes:

        task = {
            'basename' : 'Segmenting data',
            'name'     : '{0}bp'.format(window_size),
            'targets'  : [os.path.join(args.output_dir, 'segmentation_at_{0}bp.multibam'.format(window_size))],
            'file_dep' : [os.path.join(args.output_dir, 'coverage_at_{0}bp.multibam'.format(window_size))],
            }

        if args.fittings_dir:
            window_dir = os.path.join(args.output_dir, args.fittings_dir, '{0}bp'.format(window_size))
            task['actions'] = [ 'mkdir -pv {0}'.format(window_dir),
                                '%(python_exec)s %(segmentation_script)s -f {0} %(dependencies)s > %(targets)s'.format(window_dir)]
        else:
            task['actions'] = [ '%(python_exec)s %(segmentation_script)s %(dependencies)s > %(targets)s']

        yield task


def get_segmentation_bed(args):

    for window_size in args.window_sizes:
        for input_fastq in args.input_fastqs:

            yield {
                'basename' : 'Getting segmentation bed',
                'name'     : '{0}bp, {1}'.format(window_size, input_fastq),
                'actions'  : ['%(python_exec)s %(segmentation_bed_script)s -s %(dependencies)s -n {input_bam}  | tr "-" "\t" | bedtools merge > %(targets)s.unsorted'.format(input_bam=with_extension(input_fastq, '.rmdup.bam')),
                             'sort -k1,1 -k2,2n %(targets)s.unsorted > %(targets)s',
                             'rm %(targets)s.unsorted',],
                'targets'  : [with_extension(input_fastq, ".segmentation_{0}bp.bed".format(window_size))],
                'file_dep' : [os.path.join(args.output_dir, 'segmentation_at_{0}bp.multibam'.format(window_size))],
                }

def get_segmentation_bigbed(args):

    for window_size in args.window_sizes:
        for input_fastq in args.input_fastqs:

            yield {
                'basename' : 'Getting segmentation bigBed',
                'name'     : '{0}bp, {1}'.format(window_size, input_fastq),
                'actions'  : ['bedToBigBed %(dependencies)s %(genome_file)s %(targets)s',],
                'targets'  : [with_extension(input_fastq, ".segmentation_{0}bp.bb".format(window_size))],
                'file_dep' : [with_extension(input_fastq, ".segmentation_{0}bp.bed".format(window_size))],
                  }

        """
def get_matrix(args):
    last_name = []
    for window_size in reversed(args.window_sizes):

        task = {
            'basename' : 'Calculating matrix',
            'name'     : '{0}bp'.format(window_size),
            'actions'  : ['%(python_exec)s %(matrix_script)s -o -p %(num_process)s --hdf5-path %(targets)s %(dependencies)s'],
            'targets'  : [os.path.join(args.output_dir, 'segmentation_at_{0}bp.matrix.hdf5'.format(window_size))],
            'file_dep' : [os.path.join(args.output_dir, 'segmentation_at_{0}bp.multibam'.format(window_size))],
            'task_dep' : last_name
            }

        last_name = [ '{basename}:{name}'.format(**task) ]

        yield task


    return tasks_to_run
            """

class MyTaskLoader(TaskLoader):
    def __init__(self, tasks_to_run, args):
        self.tasks_to_run = tasks_to_run
        self.doit_db = os.path.join(self.get_db_dir(args), '.doit.mapping_db')
        self.num_process = args.num_process
        self.substitute_tasks(args)
        super(TaskLoader, self).__init__()

    def get_db_dir(self, args):

        output_dir = args.output_dir
        return output_dir

    def substitute_tasks(self, args):

        def get_script(script_name): return os.path.join(os.path.dirname(__file__), script_name)

        subs_dictionary = vars(args)
        subs_dictionary.update({"python_exec"         : sys.executable,
                                "segmentation_script"    : get_script('threshold_by_reads.py'),
                                "segmentation_bed_script" : get_script('extract_segmentation.py'),
                                "matrix_script"       : get_script('gam_matrix.py')})

        for task in self.tasks_to_run:
            task.options = subs_dictionary

    def load_tasks(self, cmd, opt_values, pos_args):
        config = {'verbosity': 2,
                  'dep_file' : self.doit_db,
                  'num_process' : self.num_process,
                 }
        return self.tasks_to_run, config

def main(args):

    tasks_to_run   = generate_tasks('task_map', map_fastq(args.input_fastqs))
    tasks_to_run  += generate_tasks('task_bam', sort_bam(args.input_fastqs))
    tasks_to_run  += generate_tasks('task_rmdup', remove_duplicates(args.input_fastqs))
    tasks_to_run  += generate_tasks('task_index', index_bam(args.input_fastqs))
    tasks_to_run  += generate_tasks('task_bedgraph', make_bedgraph(args.input_fastqs))
    tasks_to_run  += generate_tasks('task_bigwig', make_bigwig(args.input_fastqs))
    tasks_to_run  += generate_tasks('task_coverage', get_coverage(args))
    tasks_to_run  += generate_tasks('task_segmentation', get_segmentation(args))

    if args.bigbeds:

        tasks_to_run  += generate_tasks('task_segmentation_bed', get_segmentation_bed(args))
        tasks_to_run  += generate_tasks('task_segmentation_bigBed', get_segmentation_bigbed(args))
    #tasks_to_run  += generate_tasks('task_matrix', get_matrix(args))
    
    sys.exit(DoitMain(MyTaskLoader(tasks_to_run, args)).run([]))

if __name__ == "__main__":

    args = parser.parse_args()

    if not args.output_dir:
        args.output_dir = os.getcwd()

    main(args)
