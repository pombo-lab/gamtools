#! /usr/bin/env python

from doit.task import dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
import argparse
import os
import sys


parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-b','--input_bamfiles', metavar='INPUT_BAMFILE', required=True, nargs='+', help='One or more input bam files.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', required=True, help='File containing chromosome names and lengths')
parser.add_argument('-p','--num_process', metavar='NUM_PROCESSES', type=int, default=1, help='Number of tasks to execute in parallel')
parser.add_argument('-w','--window_sizes', metavar='WINDOW_SIZE', default=[1000,5000,10000,50000,100000,500000], type=int, nargs='+', help='File containing chromosome names and lengths')

args = parser.parse_args()

def get_coverage_tasks(input_bamfiles, window):
    """Take the path to one bamfile and return the doit tasks which would convert it to a BigWig"""

    output_dir = os.path.dirname(os.path.commonprefix(args.input_bamfiles))
    
    tasks_to_run = []

    tasks_to_run.append({
            'name' : 'Getting coverage over {0}bp windows'.format(window),
            'actions' : ['echo "window %(dependencies)s" > %(targets)s',
                         'bash -i -c "multiBamCov -bams %(dependencies)s -bed <(bedtools makewindows -w {0} -g %(genome_file)s) >> %(targets)s"' .format(window)],
            'targets' : [os.path.join(output_dir, 'coverage_at_{0}bp.multibam'.format(window))],
            'file_dep' : input_bamfiles,
           })

    return tasks_to_run


def task_to_get_segmentation(input_file, window, tasks_to_run):

    sub_vars = { 'python_exec'   : sys.executable,
             'segmentation_script' : os.path.join(os.path.dirname(__file__),'get_threshold.py'),
             'coverage_file'   : get_target(tasks_to_run, 'getting_coverage_at_{0}bp'.format(window), input_file)[0],
        }

    return {
            'name' : get_name('segmenting_at_{0}bp'.format(window), input_file),
            'actions' : [ '%(python_exec)s %(segmentation_script)s -t %(dependencies)s %(coverage_file)s > %(targets)s' % add_subs(sub_vars) ],
            'targets' : [ substitue_ext(input_file, ".segmentation_{0}bp".format(window)) ],
            'file_dep' : get_target(tasks_to_run, 'getting_saturation_at_{0}bp'.format(window), input_file),
           }

class MyTaskLoader(TaskLoader):
    def __init__(self, tasks_to_run, args):
        self.tasks_to_run = tasks_to_run
        self.doit_db = os.path.join(self.get_db_dir(args), '.doit.mapping_db')
        self.num_process = args.num_process
        self.substitute_tasks(args)
        super(TaskLoader, self).__init__()

    def get_db_dir(self, args):

        common_directory = os.path.dirname(os.path.commonprefix(args.input_bamfiles))
        return common_directory

    def substitute_tasks(self, args):

        subs_dictionary = vars(args)
        subs_dictionary.update({'dependencies' : '%(dependencies)s',
                                "targets" : '%(targets)s'})

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
    for window in args.window_sizes:
        tasks_to_run.extend( get_coverage_tasks(args.input_bamfiles, window))

    print tasks_to_run
    sys.exit(DoitMain(MyTaskLoader(tasks_to_run, args)).run([]))

if __name__ == "__main__":

    args = parser.parse_args()

    main(args)
