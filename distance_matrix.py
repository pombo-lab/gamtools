from GamTools.segmentation import map_freqs
from GamTools.cosegregation import D, Dprime, corr
import numpy as np
import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-n', '--npz_frequencies_file', required=True, help='An npz file containing co-segregation frequencies to convert to correlations')
parser.add_argument('-m', '--method', default='Dprime', help='An npz file containing co-segregation frequencies to convert to correlations')

args = parser.parse_args()

methods = { 'D' : D,
            'Dprime' : Dprime,
            'corr' : corr,
          }

correlation_file = args.npz_frequencies_file.split('.')
correlation_file[correlation_file.index('freqs')] = args.method
correlation_file = '.'.join(correlation_file)

if os.path.exists(correlation_file):
    sys.exit('File exists already, not overwriting')

method = methods[args.method]

npz_data = np.load(args.npz_frequencies_file)

windows = npz_data['windows']
freqs = npz_data['freqs']

distances = map_freqs(method, freqs)

np.savez_compressed(correlation_file, 
                    windows=windows,
                    scores=distances)
