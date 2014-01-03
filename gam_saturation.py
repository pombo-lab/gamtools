#! /usr/bin/env python

try:
    import argparse
except ImportError:
    import backport_argparse as argparse
import os
import sys
from doit.task import dict_to_task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-b','--input_bamfiles', metavar='INPUT_BAMFILE', required=True, nargs='+', help='One or more input bam files.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', required=True, help='File containing chromosome names and lengths')
parser.add_argument('-w','--window_sizes', metavar='WINDOW_SIZE', default=[1000,5000,10000,50000,100000,500000], type=int, nargs='+', help='File containing chromosome names and lengths')
parser.add_argument('-r','--resample_percentages', metavar='RESAMPLE_PERCENTAGE', default=[], type=float, nargs='+', help='List of percentages to use for resampling the bam file')

args = parser.parse_args()

def substitue_ext(target,new_ext):
    return os.path.splitext(target)[0] + new_ext

def add_subs(dictionary):
    dictionary.update({'dependencies':'%(dependencies)s',
                       'targets':'%(targets)s'})
    return dictionary

def get_name(task_name, input_file):
    return '{0}_for_{1}'.format(task_name, os.path.basename(input_file))

def get_target(task_dict, task_name, input_file):
    return task_dict[ get_name(task_name, input_file) ]['targets']

def task_to_index_bam(bam_path):
    return {
            'name' : get_name('indexing_bam', bam_path),
            'actions' : ['samtools index %(dependencies)s',],
            'targets' : [substitue_ext(bam_path,".bam.bai")],
            'file_dep' : [bam_path],
           }

def task_to_get_bam_coverage(input_file,window,resample_percentages,tasks_to_run):
    bamfiles = [ get_target(tasks_to_run, 'sorting_bam', input_file)[0] ]
    for pct in resample_percentages:

        bamfiles.append(get_target(tasks_to_run, 'resampling_bam_at_{0}pct'.format(pct), input_file)[0])

    return {
            'name' : get_name('getting_coverage_at_{0}bp'.format(window), input_file),
            'actions' : ['bash -i -c "multiBamCov -bams %(dependencies)s -bed <(bedtools makewindows -w ' + str(window) + ' -g ' + args.genome_file + ') > %(targets)s"'],
            'targets' : [substitue_ext(input_file,".coverage_{0}bp".format(window))],
            'file_dep' : bamfiles,
            'task_dep' : [get_name( 'indexing_bam', bamfiles[0])],
           }

def task_to_get_saturation(input_file,window,tasks_to_run):

    sub_vars = { 'python_exec' : sys.executable,
             'saturation_script' : os.path.join(os.path.dirname(__file__),'threshold_by_reads.py'),
             'fittings_folder' : os.path.splitext(input_file)[0] + '_%ibp_fittings' % window,
        }

    sub_vars['curve_path'] = os.path.join(sub_vars['fittings_folder'], os.path.splitext(input_file)[0] + '_saturation_' + str(window) +'.png')

    return {
            'name' : get_name('getting_saturation_at_{0}bp'.format(window), input_file),
            'actions' : [ '%(python_exec)s %(saturation_script)s -f %(fittings_folder)s -s %(curve_path)s %(dependencies)s > %(targets)s' % add_subs(sub_vars) ],
            'targets' : [substitue_ext(input_file,".saturation_{0}bp".format(window))],
            'file_dep' : get_target(tasks_to_run, 'getting_coverage_at_{0}bp'.format(window), input_file),
           }

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

def tasks_for_file(input_file):

    tasks_to_run = {}

    tasks_to_run[get_name('sorting_bam',input_file)] = {
                'name' : get_name('sorting_bam',input_file),
                'actions' : ['samtools sort %(dependencies)s %(targets)s',
                          'mv %(targets)s.bam %(targets)s'],
                'targets' : [substitue_ext(input_file,".sorted.bam")],
                'file_dep' : [input_file],
                }


    index_base_task = task_to_index_bam(get_target(tasks_to_run, 'sorting_bam', input_file)[0])
    tasks_to_run[index_base_task['name']] = index_base_task

    for pct in args.resample_percentages:
        tasks_to_run[get_name('resampling_bam_at_{0}pct'.format(pct), input_file)] = {
                    'name' : get_name('resampling_bam_at_{0}pct'.format(pct), input_file),
                    'actions' : ['samtools view -b -s ' + str(pct) + ' %(dependencies)s > %(targets)s',],
                    'targets' : [substitue_ext(input_file,".resample_{0}.bam".format(pct))],
                    'file_dep' : get_target(tasks_to_run, 'sorting_bam', input_file),
                    }

        index_resampled_task = task_to_index_bam(get_target(tasks_to_run, 'resampling_bam_at_{0}pct'.format(pct), input_file)[0])
        tasks_to_run[index_resampled_task['name']] = index_resampled_task

    for window in args.window_sizes:

        coverage_base_task = task_to_get_bam_coverage(input_file,window,args.resample_percentages,tasks_to_run)
        tasks_to_run[coverage_base_task['name']] = coverage_base_task

        saturation_base_task = task_to_get_saturation(input_file,window,tasks_to_run)
        tasks_to_run[saturation_base_task['name']] = saturation_base_task

        segmentation_base_task = task_to_get_segmentation(input_file, window, tasks_to_run)
        tasks_to_run[segmentation_base_task['name']] = segmentation_base_task

    return tasks_to_run

tasks_to_run = {}

for input_bam in args.input_bamfiles:
    tasks_to_run.update( tasks_for_file(input_bam) )

class mytaskloader(TaskLoader):
    @staticmethod
    def load_tasks(cmd, opt_values, pos_args):
        task_list = [ dict_to_task(my_task) for my_task in tasks_to_run.values() ]
        config = {'verbosity': 2}
        return task_list, config

if __name__ == "__main__":
    sys.exit(DoitMain(mytaskloader()).run([]))
