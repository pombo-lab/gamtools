from __future__ import print_function
from scipy.stats import scoreatpercentile, nbinom, norm
from scipy.optimize import fmin
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import os
import sys
import CurveFitting as cf
from bisect import bisect_left
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='plot coverage over different window sizes for a list of bam files.')
parser.add_argument('input_file', metavar='INPUT_FILE', help='Input coverage file')
parser.add_argument('-f','--fitting-folder', metavar='FITTING_FOLDER', default=False, help='If specified, save the individual curve fittings to this folder')

def cumulative_normal(x, loc, scale):
    return norm.cdf(x, loc, scale)

def cumulative_neg_binom(x, n , p):
    
    # x is in log, so transform it back
    x = list(10 ** x[1:])
    
    # Add the point 0.0
    x = [0.0] + x
        
    return nbinom.cdf(x, n, p)
    
def un_cumulative(x):
    last = None
    new_x = []
    for i in x:
        if last is None:
            last = i
            continue
        new_x.append(i - last)
        last = i
    return np.array(new_x)

def sum_to_1(x):
    return x / sum(x)

def neg_binomial(x, n, p):
    
    bin_y = cumulative_neg_binom(x, n, p)
    bin_y = un_cumulative(bin_y)
    return sum_to_1(bin_y)

def normal(x, loc, scale):

    norm_y = cumulative_normal(x, loc, scale)
    norm_y = un_cumulative(norm_y)
    return sum_to_1(norm_y)

def n_binom_plus_log_normal(params, x):
    
    bin_n, bin_p, nm_delta, nm_scale, size = params
    
    # We don't want any solutions where the mode of the lognormal is
    # less than the mode of the negative binomial. Therefore, the 
    # function is parameterized such that the position of the lognormal
    # is given as a distance from the mean of the negative binomial,
    # called nm_delta. nm_delta is always treated as positive
    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)
    
    nm_loc = np.log10(bin_mean) + np.abs(nm_delta)
    
    bin_y = neg_binomial(x, bin_n, bin_p)
    
    norm_y = normal(x, nm_loc, nm_scale)
    
    sum_y = (bin_y * (1.-abs(size))) + (norm_y * abs(size))
    
    return sum_y / sum(sum_y)

def get_fdr_threshold(x,fdr,threshold):
    threshold_index = bisect_left(fdr[::-1], threshold)
    return x[::-1][threshold_index - 1]

def mask_x_by_z(x, z):
    return [ xi for i, xi in enumerate(x) if z[i] != 0 ]

def filter_data(x, percentile, no_zeros=True):
    percentile_score = scoreatpercentile(x, percentile)
    less_than_percentile = x < percentile_score
    
    if no_zeros:
        not_a_zero = x > 0
        
        # only keep points which are both less than percentile AND not a zero
        points_to_keep = map(all,zip(less_than_percentile,not_a_zero))
        
    else:
        points_to_keep = less_than_percentile

    out_data = x[points_to_keep]

    if len(out_data):
        
        return out_data

    else:

        return x[not_a_zero]

def threshold_n_binom(params, p_value, thresh_range=range(500)):
    
    bin_n, bin_p, nm_delta, nm_scale, size = params
    
    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)
    
    nm_loc = np.log10(bin_mean) + np.abs(nm_delta)

    
    cumulative_dist = nbinom.cdf(thresh_range, bin_n, bin_p)
    
    prob_dist = sum_to_1(un_cumulative(cumulative_dist))
    index = bisect_left(prob_dist[::-1], p_value)
    return thresh_range[::-1][index]

def fit_histogram(breaks, counts, initial_params=(6.89799811e-01,   5.08503431e-01,
                                      2.69945316,  0.23982432,
                                      0.15)):
    old_stdout = sys.stdout
    with open(os.devnull, "w") as sys.stdout:
        params = fmin(cf.squared_difference, 
                                 initial_params,
                                 (n_binom_plus_log_normal, breaks, counts / float(sum(counts))),
                                 xtol=0.00001, ftol=0.00001 )
    sys.stdout = old_stdout
    return params

def get_fit_x(breaks, counts):
    
    step = (breaks[1] - breaks[0]) / 2.
    fit_x = mask_x_by_z(breaks[:-1]+step, counts)
    return fit_x

def plot_fit(breaks, counts, params):
    
    fit = sum(counts)*n_binom_plus_log_normal(params,
                                 breaks)
    
    fit_y = mask_x_by_z(fit, counts)
    
    fit_x = get_fit_x(breaks, counts)
    
    return plt.plot(fit_x, fit_y,'ro')

def plot_lognormal(breaks, counts, params):

    bin_n, bin_p, nm_delta, nm_scale, size = params
    
    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)
    
    nm_loc = np.log10(bin_mean) + np.abs(nm_delta)

    gauss_y = normal(breaks, nm_loc, nm_scale)
    gauss_y = gauss_y * sum(counts) * abs(size)
    gauss_y = mask_x_by_z(gauss_y, counts)
    
    fit_x = get_fit_x(breaks, counts)
    
    return plt.plot(fit_x, gauss_y)    
    
def plot_binom(breaks, counts, params):
    
    bin_n, bin_p, nm_delta, nm_scale, size = params
    
    bin_mean, bin_var = nbinom.stats(bin_n, bin_p)
    
    nm_loc = np.log10(bin_mean) + np.abs(nm_delta)
    
    binom_y = neg_binomial(breaks, bin_n, bin_p)
    binom_y = binom_y * sum(counts) * (1. - abs(size))
    binom_y = mask_x_by_z(binom_y, counts)
    
    fit_x = get_fit_x(breaks, counts)
    
    return plt.plot(fit_x, binom_y)
    
def plot_legend(hist_patches,
                fit_patches,
                normal_patches,
                binom_patches,
                thresh_patch):
    
    patches = (hist_patches[0],
                fit_patches[0],
                normal_patches[0],
                binom_patches[0],
                thresh_patch)
    
    labels = ('Coverage histogram',
              'Fitted values',
              'Signal distribution',
              'Noise distribution',
              'Threshold')
    
    plt.legend(patches, labels)
    
def prettify_plot(sample_name, fig):
    
    fig.suptitle('Combined fit for {0}'.format(sample_name), fontsize=18)
    
    plt.ylabel('No of windows')
    plt.xlabel('No of reads')
    locs, labels = plt.xticks()
    labels = map(int,10 ** locs)
    plt.xticks( locs, labels )    

def get_threshold(data, i, sample_name, plot_path):
    
    fig = plt.figure(figsize=( 16, 9 ))
    
    counts, breaks, hist_patches = plt.hist(np.log10(filter_data(data.iloc[:,i],99.99)),bins=50)
    
    params = fit_histogram(breaks, counts)
    
    read_threshold = threshold_n_binom(params, 0.001)

    if plot_path:

        fit_patches = plot_fit(breaks, counts, params)
        
        normal_patches = plot_lognormal(breaks, counts, params)
        
        binom_patches = plot_binom(breaks, counts, params)
        
        thresh_patch = plt.axvline(np.log10(read_threshold),color='black')
        
        plot_legend(hist_patches,
                    fit_patches,
                    normal_patches,
                    binom_patches,
                    thresh_patch)
        

        
        prettify_plot(sample_name, fig)
        
        plt.savefig(plot_path)

    plt.close()
    
    return read_threshold


if __name__ == '__main__':

    args = parser.parse_args()

    data = pd.read_csv(args.input_file,
                       delim_whitespace=True, index_col=[0,1,2])

    for i in range(len(data.columns)):

        sample_name = os.path.basename(data.columns[i]).split('.')[0]

        if args.fitting_folder:
            plot_path = os.path.join(args.fitting_folder, '{0}_fit.png'.format(sample_name))
        else:
            plot_path = None

        try:
            read_threshold = get_threshold(data, i, sample_name, plot_path)
        except IndexError:
            import pdb; pdb.set_trace()
        
        above_threshold = data.iloc[:,i] > read_threshold
        
        data.iloc[:,i] = above_threshold.astype(int)

        print('{0} done (number {1}) threshold {2}'.format(sample_name,i+1,read_threshold), file=sys.stderr)

    data.to_csv(sys.stdout, index=True, sep='\t')
