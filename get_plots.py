#! /usr/bin/env python

try:
    import argparse
except ImportError:
    import backport_argparse as argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-i','--input_segmentations', metavar='INPUT_SEGMENT_FILE', required=True, nargs='+', help='One or more input segmentation files.')
parser.add_argument('-o','--output_dir', metavar='OUTPUT_DIR', required=True)

args = parser.parse_args()

def read_files(input_files):
    temp_matrix = []
    for path in input_files:
        col_name = os.path.basename(path).split('_')[0]
        temp_matrix.append(pd.read_csv(path,delim_whitespace=True,names=[col_name],index_col=[0,1]))
        
    matrix = pd.concat(temp_matrix,axis=1)
    return matrix

def get_num(df,index1,index2,val1,val2):
    return len(df[(df[index1] == val1) & (df[index2] == val2)])

def get_table(df, index1, index2):
    return np.array([[ get_num(df, index1, index2, False, False), get_num(df, index1, index2, False, True)],
                    [ get_num(df, index1, index2, True, False), get_num(df, index1, index2, True, True)]]) + 0.01

def odds_ratio(N):
    return np.log(N[0,0]) + np.log(N[1,1]) - np.log(np.log(N[0,1])) - np.log(N[1,0])

def get_chr_heatmap(df,chrom):
    data = df.ix[chrom].transpose()
    img_data = []
    total = len(data.columns) **2
    for i,col1 in enumerate(data.columns):
        img_data.append([])
        for q,col2 in enumerate(data.columns):
            n = (i + 1) * (q + 1)
            if n % (total / 10) == 0:
                print 's1','.'
            img_data[i].append(odds_ratio(get_table(data,col1,col2)))
    return np.array(img_data)

if __name__ == "__main__":

    if not os.path.isdir(args.output_dir):
        print 'dir does not exist'
        os.mkdir(args.output_dir)
    full_m = read_files(args.input_segmentations)

    names = set(map(lambda x:x[0],full_m.index))
    print 's2',names
    for chrm in names:
        if chrm[-7:] == '_random':
            continue
        print 's3','getting heatmap for %s' % chrm
        img_data = get_chr_heatmap(full_m, chrm)
        plt.figure(figsize=(20,20))
        plt.imshow(img_data, interpolation='none')
        plt.savefig(os.path.join(args.output_dir,'%s.png' % chrm))
