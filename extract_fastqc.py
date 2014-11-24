#!/usr/bin/python

# Shamelessly stolen from https://code.google.com/p/bioinformatics-misc/source/browse/trunk/fastqc_to_pgtable.py?spec=svn93&r=93

import argparse
import numpy as np
import re

parser = argparse.ArgumentParser(description= """
    Convert the FastQC report (file: fastqc_data.txt) to a single line with columns
    tab-separated. The input is the directory produced by fastqc which must
    contain the file fastqc_data.txt with the actual stats.
    """)

parser.add_argument('--infiles', '-i',
                   nargs='+',
                   required= True,
                   help='''Fastqc *directory* where to look for file fastqc_data.txt.
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------

def parse_module(fastqc_module):
    """
    Parse a fastqc module from the table format to a line format (list).
    Input is list containing the module. One list-item per line. E.g.:

    fastqc_module= [
        '>>Per base sequence quality    pass',
        '#Base  Mean    Median  Lower Quartile  Upper Quartile  10th Percentile 90th Percentile',
        '1  36.34   38.0    35.0    39.0    31.0    40.0',
        '2  35.64   38.0    35.0    39.0    28.0    40.0',
        '3  35.50   38.0    34.0    39.0    28.0    40.0',
        ...
        ]
    
    Return a list like this where each sublist after 1st is a column:
    ['pass', ['1', '2', '3', ...], ['36,34', '35.64', '35.50', ...], ['40.0', '40.0', '40.0', ...], ...]
    """
    row_list= []
    module_header= fastqc_module[0]
    module_name= module_header.split('\t')[0]
    row_list.append(module_header.split('\t')[1]) ## Line with module name >>Per base ...    pass/warn/fail    
    
    # Handle odd cases:
    # Where no table is returned:
    if len(fastqc_module) == 1 and module_name == '>>Overrepresented sequences':
        return(row_list + [[]]*4)
    if len(fastqc_module) == 1 and module_name == '>>Kmer Content':
        return(row_list + [[]]*5)
    # Table is not the secod row:
    if module_name == '>>Sequence Duplication Levels':
        tot_dupl= fastqc_module[1].split('\t')[1]
        row_list.append(tot_dupl)
        del(fastqc_module[1])
   
    ## Conevrt table to list of lists:
    tbl= []
    for line in fastqc_module[2:]:
        tbl.append(line.split('\t'))
    ## Put each column in a list:
    nrows= len(tbl)
    ncols= len(tbl[0])
    for i in range(0, ncols):
        col= []
        for j in range(0, nrows):
            col.append(tbl[j][i])
        row_list.append(col)
    return(row_list)
 
def is_mono_repeat(kmer):
    if len(set(kmer)) == 1:
        return True
    else:
        return False

def is_di_repeat(kmer):
    if not len(set(kmer)) == 2:
        return False

    if len(set(kmer[::2])) == 1 and len(set(kmer[1::2])) == 1:
        return True
    else:
        return False

def get_kmer_summary(module):

    kmer_data = parse_module(module)
    kmers = kmer_data[1]
    counts = map(float, kmer_data[3])

    summary_data = {'dinucleotide_repeats':0,
                    'mononucleotide_repeats':0}

    for kmer, count in zip(kmers, counts):
        if is_mono_repeat(kmer):
            summary_data['mononucleotide_repeats'] += count
        elif is_di_repeat(kmer):
            summary_data['dinucleotide_repeats'] += count

    return summary_data

def get_avg_qual(module):

    qual_data = parse_module(module)
    qualities, counts = np.array(map(int, qual_data[1])), np.array(map(float, qual_data[2]))
    avg = (qualities * counts).sum() / counts.sum()
    return { 'avg_quality':avg }

def get_sample(filename):

    try:
        return re.search('F[0-9]{1,2}[A-H][0-9]', filename).group(0)
    except AttributeError:
        return re.search('E[0-9]', filename).group(0)

def process_file(filename):
    fq= open(filename).readlines()
    fq= [x.strip() for x in fq]

    fastqc_dict = {}

    ## Get start and end position of all modules:
    mod_start= []
    mod_end= []
    for i in range(0, len(fq)):
        line= fq[i]

        if line == '>>END_MODULE':
            mod_end.append(i)
        elif line.startswith('>>'):
            mod_start.append(i)
        else:
            pass

    ## Start processing modules. First one (Basic statitics) is apart:
    for s, e in zip(mod_start[1:], mod_end[1:]):
        module= fq[s:e]
        module_name= module[0].split('\t')[0]
        if module_name == '>>Kmer Content':

            fastqc_dict.update(get_kmer_summary(module))

        if module_name == '>>Per sequence quality scores':
            
            fastqc_dict.update(get_avg_qual(module))

    print get_sample(filename), fastqc_dict['avg_quality'], fastqc_dict['mononucleotide_repeats'], fastqc_dict['dinucleotide_repeats']

print 'Sample   Avg_quality   mononucleotide_repeats   dinucleotide_repeats'

for filename in args.infiles:
    process_file(filename)
