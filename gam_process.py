#! /usr/bin/env python

from wrapit.api import run
from wrapit.parser import parser
import os
import sys


parser.description='Map GAM-seq reads and create BigWigs and fastQC files'
parser.add_argument('input_fastqs', metavar='INPUT_FASTQ', nargs='+', help='One or more input fastq files.')
parser.add_argument('-g','--genome_file', metavar='GENOME_FILE', required=True, 
                    help='File containing chromosome names and lengths')
parser.add_argument('-o','--output_dir', metavar='OUPUT_DIRECTORY', help='Write segmentation, matrix etc. to this directory')
parser.add_argument('-f','--fittings_dir', metavar='FITTINGS_DIR', help='Write segmentation curve fitting plots to this directory')
parser.add_argument('-i','--bigwigs', action='append_const',
                    dest='to_run', const='Converting bedgraph to bigwig', help='Make bigWig files.')
parser.add_argument('-b','--bigbeds', action='append_const',
                    dest='to_run', const='Getting segmentation bigBed', help='Make bed files for segmentation')
parser.add_argument('-c','--do-qc', action='append_const',
                    dest='to_run', const='do_qc', help='Run various qc scripts.')
parser.add_argument('-w','--window_sizes', metavar='WINDOW_SIZE', default=[1000,5000,10000,50000,100000,500000], type=int, nargs='+', help='File containing chromosome names and lengths')
parser.add_argument('--qc-window-size', type=int, help='Use this window size for qc (default is median window size).')
parser.add_argument('--additional-qc-files', nargs='*', default=[],
                    help='Any additional qc files to filter on')
parser.add_argument('-q','--minimum-mapq', metavar='MINIMUM_MAPQ', default=0, type=int, help='Filter out any mapped read with a mapping quality less than x (default is not to filter based on quality')

def get_middle_value(_list):
    return sorted(_list)[len(_list)/2]

def get_script(script_name): return os.path.join(os.path.dirname(__file__), script_name)

parser.set_defaults(python_exec=sys.executable,
                    segmentation_script=get_script('threshold_by_reads.py'),
                    segmentation_bed_script=get_script('extract_segmentation.py'),
                    matrix_script=get_script('gam_matrix.py'),
                    contamination_script=get_script('extract_contamination.py'),
                    quality_script=get_script('extract_fastqc.py'),
                    segmentation_stats_script=get_script('segmentation_stats.py'),
                    stats_merge_script=get_script('merge_stats.py'),
                    pass_qc_script=get_script('pass_qc.py'),
                    filter_script=get_script('select_samples.py'),
                    default_stats=['contamination_stats.txt', 'mapping_stats.txt',
                                   'quality_stats.txt', 'segmentation_stats.txt'],
                    to_run=['Segmenting data'])

def with_extension(input_fastq, extension): 
    output_dir = os.path.dirname(input_fastq)
    basename = os.path.splitext(os.path.basename(input_fastq))[0]
    return os.path.join(output_dir, basename + extension)

def task_map_fastq(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Mapping fastq',
                "actions"  : ['bowtie2 -x genome -U %(dependencies)s | sed \'/XS : /d\' | samtools view -q %(minimum_mapq)s -F 4 -bS - > %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".bam")],
                "file_dep" : [input_fastq],
              }

def task_sort_bam(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Sorting bam',
                "actions"  : ['samtools sort %(dependencies)s %(targets)s',
                              'mv %(targets)s.bam %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".sorted.bam")],
                "file_dep" : [with_extension(input_fastq, ".bam")],
              }

def task_remove_duplicates(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Removing duplications',
                "actions"  : ['samtools rmdup -s %(dependencies)s %(targets)s',],
                "targets"  : [with_extension(input_fastq, ".rmdup.bam")],
                "file_dep" : [with_extension(input_fastq, ".sorted.bam")],
              }

def task_index_bam(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Indexing bam',
                "actions"  : ['samtools index %(dependencies)s',],
                "targets"  : [with_extension(input_fastq, ".rmdup.bam.bai")],
                "file_dep" : [with_extension(input_fastq, ".rmdup.bam")],
              }

def task_make_bedgraph(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Converting bam to bedgraph',
                "actions"  : ['genomeCoverageBed -bg -ibam %(dependencies)s -g %(genome_file)s > %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".bedgraph")],
                "file_dep" : [with_extension(input_fastq, ".rmdup.bam")],
              }

def task_make_bigwig(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Converting bedgraph to bigwig',
                "actions"  : ['bedGraphToBigWig %(dependencies)s %(genome_file)s %(targets)s'],
                "targets"  : [with_extension(input_fastq, ".bw")],
                "file_dep" : [with_extension(input_fastq, ".bedgraph")],
              }

def task_get_coverage(args):
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

def task_get_segmentation(args):
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

def task_get_segmentation_bed(args):

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

def task_get_segmentation_bigbed(args):

    for window_size in args.window_sizes:
        for input_fastq in args.input_fastqs:

            yield {
                'basename' : 'Getting segmentation bigBed',
                'name'     : '{0}bp, {1}'.format(window_size, input_fastq),
                'actions'  : ['bedToBigBed %(dependencies)s %(genome_file)s %(targets)s',],
                'targets'  : [with_extension(input_fastq, ".segmentation_{0}bp.bb".format(window_size))],
                'file_dep' : [with_extension(input_fastq, ".segmentation_{0}bp.bed".format(window_size))],
                  }

def fastqc_data_file(input_fastq):
    
    base_folder = '.'.join(input_fastq.split('.')[:-1])

    fastqc_folder = base_folder + '_fastqc'

    return os.path.join(fastqc_folder, 'fastqc_data.txt')

def task_run_fastqc(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Running fastqc',
                "actions"  : ['fastqc %(dependencies)s',
                              'head %(targets)s'],
                "targets"  : [fastqc_data_file(input_fastq)],
                "file_dep" : [input_fastq],
              }

def task_run_fastq_screen(args):
    for input_fastq in args.input_fastqs:
        yield {
                "name"     : input_fastq,
                "basename" : 'Running fastq_screen',
                "actions"  : ['fastq_screen %(dependencies)s'],
                "targets"  : [input_fastq + '_screen.txt'],
                "file_dep" : [input_fastq],
              }

def task_extract_contamination(args):
    fastq_screen_files = [ input_fastq + '_screen.txt' for input_fastq in args.input_fastqs ]
    return {
            'actions' : ['%(python_exec)s %(contamination_script)s %(dependencies)s > %(targets)s'],
            'targets' : [os.path.join(args.output_dir, 'contamination_stats.txt')],
            'file_dep' : fastq_screen_files,
           }

def task_extract_quality(args):
    fastqc_files = [ fastqc_data_file(input_fastq) for input_fastq in args.input_fastqs ]
    return {
            'actions' : ['%(python_exec)s %(quality_script)s -i %(dependencies)s > %(targets)s'],
            'targets' : [os.path.join(args.output_dir, 'quality_stats.txt')],
            'file_dep' : fastqc_files,
           }

def task_mapping_stats(args):
    return {
            'actions' : ['echo "Sample\tReads_sequenced\tReads_mapped\tUnique_reads_mapped" > %(targets)s',
                         """for fq in %(dependencies)s; 
                         do 
                            sample=$(basename $fq);
                            nlines=$(zcat ${fq} | wc -l);
                            nreads=$(echo "${nlines} / 4" | bc);
                            nmap=$(samtools view -c ${fq%%.*}.bam);
                            nnodup=$(samtools view -c ${fq%%.*}.rmdup.bam);
                            echo "${sample%%%%.*}\t${nreads}\t${nmap}\t${nnodup}"; 
                         done >> %(targets)s"""],
            'targets' : [os.path.join(args.output_dir, 'mapping_stats.txt')],
            'file_dep' : args.input_fastqs,
           }

def task_segmentation_stats(args):

    input_segmentation = os.path.join(args.output_dir, 'segmentation_at_{0}bp.multibam'.format(args.qc_window_size))

    return {
            'actions' : ['%(python_exec)s %(segmentation_stats_script)s %(dependencies)s > %(targets)s'],
            'targets' : [os.path.join(args.output_dir, 'segmentation_stats.txt')],
            'file_dep' : [input_segmentation]
           }

def task_merge_stats(args):

    files_to_merge = [os.path.join(args.output_dir, sf) for sf in args.default_stats] + args.additional_qc_files

    return {
            'actions' : ['%(python_exec)s %(stats_merge_script)s %(dependencies)s > %(targets)s'],
            'targets' : [os.path.join(args.output_dir, 'merged_stats.txt')],
            'file_dep' : files_to_merge
           }

def task_samples_passing_qc(args):

    return {
            'actions' : ['%(python_exec)s %(pass_qc_script)s %(dependencies)s > %(targets)s'],
            'targets' : [os.path.join(args.output_dir, 'pass_qc.txt')],
            'file_dep' : [os.path.join(args.output_dir, 'merged_stats.txt')]
           }

def task_filter_segmentation(args):
    for window_size in args.window_sizes:

        segmentation_file = os.path.join(args.output_dir, 'segmentation_at_{0}bp.passqc.multibam'.format(window_size))
        task = {
            'basename' : 'Filtering out failed qc samples',
            'name'     : '{0}bp'.format(window_size),
            'targets'  : [os.path.join(args.output_dir, 'segmentation_at_{0}bp.passqc.multibam'.format(window_size))],
            'file_dep' : [os.path.join(args.output_dir, 'pass_qc.txt')],
            'actions'  : [ 'cat %(dependencies)s | xargs %(python_exec)s %(filter_script)s -s {seg_file} -g -n > %(targets)s'.format(seg_file=segmentation_file)]
            }

        yield task

def task_do_qc():
    return {'actions': None,
            'task_dep': ['Filtering out failed qc samples']}

def main(args):

    run(globals(), args, args.to_run)

if __name__ == "__main__":

    args = parser.parse_args()
    if args.qc_window_size is None:
        args.qc_window_size = get_middle_value(args.window_sizes)

    if not args.output_dir:
        args.output_dir = os.getcwd()

    main(args)
