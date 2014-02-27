import pandas as pd
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='Extract from a segmentation file only the positive windows for a particular sample')
parser.add_argument('-s1','--segmentation-file1', metavar='SEGMENTATION_FILE_1', required=True, help='A file containing the segmentation of all samples')
parser.add_argument('-s2','--segmentation-file2', metavar='SEGMENTATION_FILE_2', required=True, help='A file containing the segmentation of all samples')
parser.add_argument('-n','--sample-name', metavar='SAMPLE_NAME', required=True, help='Name of the sample to extract')

args = parser.parse_args()

data1 = pd.read_csv(args.segmentation_file1, sep='\s*')
data2 = pd.read_csv(args.segmentation_file2, sep='\s*')
data1['start'], data1['stop'] = zip(*map(lambda x: x.split('-'),data1['window']))
data1['name'] = args.sample_name
data1['score'] = 0
data1['strand'] = '+'
data1['tstart'] = data1['start']
data1['tstop'] = data1['stop']

a = data1[args.sample_name]
b = data2[args.sample_name]

seg_choices = {(0, 0):0, (0, 1):1, (1, 0):1, (1, 1):2}
col_choices = {(0, 0):'255,255,255',
               (0, 1):'255,0,0',
               (1, 0):'255,0,0',
               (1, 1):'0,0,0'}

data1['col'] = [ col_choices[x] for x in zip(a,b)]
data1['seg'] = [ seg_choices[x] for x in zip(a,b)]

subset = data1.loc[data1['seg'] > 0,['chrom','start','stop','name','score','strand','tstart','tstop','col']]

sample_name = os.path.basename(args.sample_name).split('.')[0]
print 'track name="{0}" description="Sample {0} segmentation difference" itemRgb="On"'.format(sample_name)
subset.to_csv(sys.stdout, header=False, index=False, sep='\t')
