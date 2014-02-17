import numpy as np
import pandas as pd
import os
from bisect import bisect_left
import argparse
from scipy.stats import linregress
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description='Given a coverage file and a segmentation file for a saturation curve, print some information about whether sampling has neared saturation')
parser.add_argument('-b','--input_coverage_files', metavar='INPUT_COVERAGE', required=True, nargs='+', help='One or more input coverage multibams.')

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

def get_cov_seg_data(subs_path):
    cpath = subs_path.format('coverage')
    spath = subs_path.format('segmentation')
    cdata = pd.read_csv(cpath, delim_whitespace=True, index_col=[0,1,2], header=None)
    nwindow = len(cdata.index)
    sdata = pd.read_csv(spath, delim_whitespace=True, index_col=[0,1], header=None, skiprows=1)
    csum, ssum = list(cdata.sum()), list(sdata.sum())
    csum, ssum = zip(*[ i for i in zip(csum,ssum) if i[1] != 0.0])
    csum, ssum = np.array([0.0] + list(csum)), np.array([0.0] + list(ssum))
    ssum = ssum / float(nwindow)
    
    #sort coverage and saturation data by coverage depth
    csum, ssum = zip(*sorted(zip(csum, ssum),key=lambda i: i[0]))
    return csum, ssum

def get_curve_split(csum, ssum, s):
    whole = linregress(csum, ssum)[2]
    fps = linregress(csum[:s], ssum[:s])
    fp = fps[2]
    sps = linregress(csum[s:], ssum[s:])
    sp = sps[2]
    mean_halves = np.mean([fp,sp])
    if fps[0] > sps[0]:
        return (mean_halves - whole) * 10
    else:
        return 0.0

def get_curve_index(csum,ssum):
    curve_scores = []
    for s in range(2, len(csum) - 1):
        curve_scores.append(get_curve_split(csum, ssum, s))
    return max(curve_scores)

def get_saturation_point(csum, ssum):
    sat_point = 0.85 * max(ssum)
    sat_ix = bisect_left(ssum, sat_point)
    return csum[sat_ix], ssum[sat_ix]
    
def plot_saturation(csum, ssum, subs_path):
    plt.plot(csum, ssum, 'bo')
    plt.title(os.path.basename(subs_path).split('.')[0])

    curve_index = get_curve_index(csum, ssum)
    
    if curve_index > 0.4:
        sat_x, sat_y = get_saturation_point(csum, ssum)
        plt.axvline(sat_x, color='black')
        return sat_x, sat_y, curve_index
    else:
        return 'NA', 'NA', curve_index
        
def get_saturation(subs_path, img_path):
    
    csum, ssum = get_cov_seg_data(subs_path)

    sat_x, sat_y, curve_index = plot_saturation(csum, ssum, subs_path)
    
    if sat_x == 'NA':
        saturated = False
    else:
        saturated = True
    
    print os.path.basename(subs_path).split('.')[0], saturated, curve_index, sat_x, sat_y  

    plt.savefig(img_path)
    plt.cla()

def main(input_coverage_files):
    for f in input_coverage_files:
        f = f.replace('coverage', '{0}')
        img_path = f.format('saturation') + '.png'
        get_saturation(f, img_path)

if __name__ == "__main__":

    args = parser.parse_args()

    main(args.input_coverage_files)
