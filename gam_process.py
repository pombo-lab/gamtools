#! /usr/bin/env python

from doit.task import dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
import argparse
import os
import sys


parser = argparse.ArgumentParser(description='Map GAM-seq reads and create BigWigs and fastQC files')
parser.add_argument('input_fastqs', metavar='INPUT_FASTQ', nargs='*', help='One or more input fastq files.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', help='File containing chromosome names and lengths')
parser.add_argument('-p','--num_process', metavar='NUM_PROCESSES', type=int, default=1, help='Number of tasks to execute in parallel')
parser.add_argument('-o','--output_dir', metavar='OUPUT_DIRECTORY', help='Write segmentation, matrix etc. to this directory')
parser.add_argument('-w','--window_sizes', metavar='WINDOW_SIZE', default=[1000,5000,10000,50000,100000,500000], type=int, nargs='+', help='File containing chromosome names and lengths')

args = parser.parse_args()

def get_mapping_tasks(input_fastq):
    """Take the path to one bamfile and return the doit tasks which would convert it to a BigWig"""

    output_dir = os.path.dirname(input_fastq)
    basename = os.path.splitext(os.path.basename(input_fastq))[0]

    def with_extension(extension): return os.path.join(output_dir, basename + extension)

    tasks_to_run = []

    tasks_to_run.append({
                "name"     : 'Mapping fastq',
                "actions"  : ['bowtie2 -x genome -U %(dependencies)s | sed \'/XS : /d\' | samtools view -F 4 -bS - > %(targets)s'],
                "targets"  : [with_extension(".bam")],
                "file_dep" : [input_fastq],
                })

    tasks_to_run.append({
                "name"     : 'Sorting bam',
                "actions"  : ['samtools sort %(dependencies)s %(targets)s',
                              'mv %(targets)s.bam %(targets)s'],
                "targets"  : [with_extension(".sorted.bam")],
                "file_dep" : [with_extension(".bam")],
                })

    tasks_to_run.append({
                "name"     : 'Removing duplications',
                "actions"  : ['samtools rmdup -s %(dependencies)s %(targets)s',],
                "targets"  : [with_extension(".rmdup.bam")],
                "file_dep" : [with_extension(".sorted.bam")],
                })

    tasks_to_run.append({
                "name"     : 'Converting bam to bedgraph',
                "actions"  : ['genomeCoverageBed -bg -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
                "targets"  : [with_extension(".bedgraph")],
                "file_dep" : [with_extension(".rmdup.bam")],
                })

    tasks_to_run.append({
                "name"     : 'Converting bedgraph to bigwig',
                "actions"  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                "targets"  : [with_extension(".bw")],
                "file_dep" : [with_extension(".bedgraph")],
                })

    tasks_to_run.append({
                'name'     : 'Indexing bam',
                'actions'  : ['samtools index %(dependencies)s',],
                'targets'  : [with_extension(".rmdup.bam.bai")],
                'file_dep' : [with_extension(".rmdup.bam")],
                })

    for task in tasks_to_run:
        task["name"] = task["name"] + " (%s)" % task["targets"][0]

    return tasks_to_run

def get_coverage_tasks(tasks_to_run, window, output_dir):
    """Take the path to one bamfile and return the doit tasks which would convert it to a BigWig"""

    bamfiles = []
    for task in tasks_to_run:
        if task['name'].startswith('Removing'):
            bamfiles.extend(task['targets'])
    
    indexes = []
    for task in tasks_to_run:
        if task['name'].startswith('Indexing'):
            indexes.append(task['name'])

    tasks_to_run.append({
            'name' : 'Getting coverage over {0}bp windows'.format(window),
            'actions' : ['echo "window %(dependencies)s" > %(targets)s',
                         'bash -i -c "multiBamCov -bams %(dependencies)s -bed <(bedtools makewindows -w {0} -g %(genome_file)s) >> %(targets)s"' .format(window)],
            'targets' : [os.path.join(output_dir, 'coverage_at_{0}bp.multibam'.format(window))],
            'file_dep' : bamfiles,
            'task_dep' : indexes,
           })

    tasks_to_run.append({
            'name'     :  'Getting threshold at {0}bp'.format(window),
            'actions'  : ['%(python_exec)s %(threshold_script)s --header-lines 1 %(dependencies)s > %(targets)s'],
            'targets'  : [os.path.join(output_dir, 'threshold_at_{0}bp.txt'.format(window))],
            'file_dep' : [os.path.join(output_dir, 'coverage_at_{0}bp.multibam'.format(window))],
            })

    tasks_to_run.append({
            'name'     : 'Segmenting data at {0}bp'.format(window),
            'actions'  : ['%(python_exec)s %(segmentation_script)s --header-lines 1 -t %(dependencies)s {coverage_file} > %(targets)s'.format(
                                                         coverage_file=os.path.join(output_dir, 'coverage_at_{0}bp.multibam'.format(window)))],
            'targets'  : [os.path.join(output_dir, 'segmentation_at_{0}bp.multibam'.format(window))],
            'file_dep' : [os.path.join(output_dir, 'threshold_at_{0}bp.txt'.format(window))],
            })

    """
    tasks_to_run.append({
            'name'     :  'Getting matrix at {0}bp'.format(window),
            'actions'  : ['%(python_exec)s %(matrix_script)s --hdf5-path %(targets)s %(dependencies)s'],
            'targets'  : [os.path.join(output_dir, '{0}bp_matrix.hdf5'.format(window))],
            'file_dep' : [os.path.join(output_dir, 'segmentation_at_{0}bp.multibam'.format(window))],
            })
            """

    return tasks_to_run

def get_segmentation_bed_tasks(tasks_to_run, segmentation_file, window):

    bamfiles = []
    for task in tasks_to_run:
        if task['name'].startswith('Removing'):
            bamfiles.extend(task['targets'])

    for input_bam in bamfiles:

        output_dir = os.path.dirname(input_bam)
        basename = os.path.splitext(os.path.basename(input_bam))[0]

        def with_extension(extension): return os.path.join(output_dir, basename + extension)

        tasks_to_run.append({
                    'name'     : 'Getting segmentation bed at {0}bp for {1}'.format(window, input_bam),
                    'actions'  : ['%(python_exec)s %(segmentation_bed_script)s -s %(dependencies)s -n {input_bam}  | tr "-" "\t" | bedtools merge > %(targets)s'.format(input_bam=input_bam),],
                    'targets'  : [with_extension(".segmentation_{0}bp.bed".format(window))],
                    'file_dep' : [segmentation_file],
                    })

        tasks_to_run.append({
                    'name'     : 'Sorting segmentation bed at {0}bp for {1}'.format(window, input_bam),
                    'actions'  : ['sort -k1,1 -k2,2n %(dependencies)s > %(targets)s',],
                    'targets'  : [with_extension(".segmentation_{0}bp.sorted.bed".format(window))],
                    'file_dep' : [with_extension(".segmentation_{0}bp.bed".format(window))],
                    })

        tasks_to_run.append({
                    'name'     : 'Making segmentation bigBed at {0}bp for {1}'.format(window, input_bam),
                    'actions'  : ['bedToBigBed %(dependencies)s %(genome_file)s %(targets)s',],
                    'targets'  : [with_extension(".segmentation_{0}bp.bb".format(window))],
                    'file_dep' : [with_extension(".segmentation_{0}bp.sorted.bed".format(window))],
                    })

    return tasks_to_run

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
        subs_dictionary.update({'dependencies'        : '%(dependencies)s',
                                "targets"             : '%(targets)s',
                                "python_exec"         : sys.executable,
                                "threshold_script"    : get_script('threshold_by_reads.py'),
                                "segmentation_script" : get_script('get_threshold.py'),
                                "segmentation_bed_script" : get_script('extract_segmentation.py'),
                                "matrix_script"       : get_script('gam_matrix.py')})

        for task in self.tasks_to_run:
            task['actions'] = [ action % subs_dictionary for action in task['actions'] ]

    def load_tasks(self, cmd, opt_values, pos_args):
        task_list = [ dict_to_task(my_task) for my_task in self.tasks_to_run ]
        config = {'verbosity': 2,
                  'dep_file' : self.doit_db,
                  'num_process' : self.num_process,
                 }
        return task_list, config

def main(args):

    tasks_to_run = []
    for input_fastq in args.input_fastqs:
        tasks_to_run.extend( get_mapping_tasks(input_fastq))

    for window in args.window_sizes:
        tasks_to_run = get_coverage_tasks(tasks_to_run, window, args.output_dir)
        segmentation_file = os.path.join(args.output_dir, 'segmentation_at_{0}bp.multibam'.format(window))

        tasks_to_run = get_segmentation_bed_tasks(tasks_to_run, segmentation_file, window)
    
    sys.exit(DoitMain(MyTaskLoader(tasks_to_run, args)).run([]))

if __name__ == "__main__":

    args = parser.parse_args()

    if not args.output_dir:
        args.output_dir = os.getcwd()

    main(args)
