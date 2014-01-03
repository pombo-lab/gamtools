from scipy.stats import scoreatpercentile
from scipy.optimize import fmin
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
import CurveFitting as cf
from CurveFitting import NoMaximaException
from bisect import bisect_left
try:
    import argparse
except ImportError:
    import backport_argparse as argparse
import itertools

parser = argparse.ArgumentParser(description='plot coverage over different window sizes for a list of bam files.')
parser.add_argument('input_file', metavar='INPUT_FILE', help='Input coverage file')
parser.add_argument('-f','--fitting-folder', metavar='FITTING_FOLDER', default=False, help='If specified, save the individual curve fittings to this folder')
parser.add_argument('-s','--saturation-curve', metavar='CURVE_FILENAME', help='If specified, save the saturation curve to this file')
parser.add_argument('--header-lines', metavar='NUMBER_OF_LINES', default=0, type=int, help='Skip this number of lines in the header')
args = parser.parse_args()

def get_kde(data):
    kernel = cf.gaussian_kde(data)
    kernel.silverman_factor()
    data_max = max(data)
    step = data_max / 300
    points = np.arange(0,data_max,step)
    return points,kernel.evaluate(points)

def get_data(fobj):
    data = []
    for line in itertools.islice(fobj, args.header_lines, None):
        fields = line.split()
        data.append(map(int,fields[3:]))
    return np.transpose(np.array(data))

def get_fdr_threshold(x,fdr,threshold):
    threshold_index = bisect_left(fdr[::-1], threshold)
    return x[::-1][threshold_index - 1]

def do_fitting(data,stringency=5,filter_pct=99.9):
    filtered_data = data[data>0]
    filtered_data = filtered_data[filtered_data<scoreatpercentile(filtered_data,filter_pct)]
    filtered_data = np.log10(filtered_data)
    kde_points = get_kde(filtered_data)
    kde_x, kde_y = kde_points
    try:
        last_max = cf.recursive_last_maximum(kde_x, kde_y, stringency)
    except NoMaximaException:
        print 'Couldnt find last_max'
        last_max = len(kde_x) - 1
        return None, filtered_data, kde_points, last_max
    unfit_segment, fit_segment = cf.split_line(kde_x, kde_y, last_max)
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")
    fit_results = fmin(cf.squared_difference, np.array([ kde_x[last_max],  0.38570735,  kde_y[last_max]]),(cf.half_gaussian,) + fit_segment)
    sys.stdout = old_stdout
    
    return fit_results, filtered_data, kde_points, last_max

def get_threshold(fit_results,kde_points,fdr_thresh=0.10):
    
    kde_x, kde_y = kde_points
    z = cf.full_gaussian(np.array(fit_results),kde_x)
    fdr = cf.calculate_fdr(z, kde_y)
    fdr_thresh = get_fdr_threshold(kde_x,fdr,fdr_thresh)
    return fdr_thresh

def get_coverage(unfiltered_data, threshold):
    return len(unfiltered_data[unfiltered_data>threshold]) / float(len(unfiltered_data))

def plot_both(fname,filtered_data,fit_results,kde_points,fdr_thresh,last_max,title=False):
    
    fig = plt.figure(figsize=( 16, 9 ))
    ax1 = fig.add_subplot(111)
    n,bins,hist_patches = ax1.hist(filtered_data,100)
    max_n = max(n)
    
    kde_x, kde_y = kde_points
    scaling_factor = ( max_n / max(kde_y))
    kde_y = kde_y * scaling_factor
    unfit_segment, fit_segment = cf.split_line(kde_x, kde_y, last_max)
    
    fit_patch = ax1.plot(*fit_segment,color='r')
    unfit_patch = ax1.plot(*unfit_segment,color='g')
    
    z = cf.full_gaussian(np.array(fit_results),kde_x)
    z = z * scaling_factor
    result_gaussian_patch = ax1.plot(kde_x,z,'purple')  
    
    thresh_patch = ax1.axvline(fdr_thresh,color='black')
    if title:
         fig.suptitle(title, fontsize=18)
    ax1.set_ylabel('No of windows')
    ax1.set_xlabel('No of reads')
    ax1.set_xticklabels( 10 ** ax1.get_xticks() )
    
    plt.figlegend((hist_patches[0],fit_patch[0],unfit_patch[0],result_gaussian_patch[0],thresh_patch),
    ('Depth Histogram','KDE (fitted)','KDE (unfit)','Fitting Result','Threshold'),loc = 'upper right')
    fig.savefig(fname)

def plot_saturation(fname,counts,coverages,title=False):
    fig = plt.figure(figsize=( 12, 7 ))
    ax1 = fig.add_subplot(111)
    plt.plot(counts,coverages,'ro')
    ax1.set_ylabel('% of windows over threshold')
    ax1.set_xlabel('% of original reads')
    if title:
         fig.suptitle(title, fontsize=18)
    ax1.set_xlim(left=0.0)
    ax1.set_ylim(bottom=0.0)
    fig.savefig(fname)

if args.fitting_folder and not os.path.isdir(args.fitting_folder):
    os.mkdir(args.fitting_folder)

data = get_data(open(args.input_file))

counts = np.array(map(sum,data))
max_count = float(max(counts))
counts = counts / max_count

coverages = []
for i,count in enumerate(counts):
    
    fit_results, filtered_data, kde_points, last_max = do_fitting(data[i])
    kde_x, kde_y = kde_points
    
    if not fit_results is None:
        fdr_thresh = get_threshold(fit_results,kde_points)
    else:
        fdr_thresh = 1.5
    if fdr_thresh < np.log10(3): fdr_thresh = np.log10(3)
    coverage = get_coverage(data[i], 10**fdr_thresh)
    coverages.append(coverage)
    print sum(data[i]), coverage, 10**fdr_thresh
    if args.fitting_folder:
        filename = '%s_resample_%i%%.png' % (os.path.basename(args.input_file), 100 * count)
        title = '%s resampled at %i%%.png' % (os.path.basename(args.input_file), 100 * count)
        plot_both(os.path.join(args.fitting_folder, filename), filtered_data,fit_results,kde_points,fdr_thresh,last_max,title)

if args.saturation_curve:
    plot_saturation(args.saturation_curve, counts, coverages)
