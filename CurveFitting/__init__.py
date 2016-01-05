import numpy as np
from scipy.stats import gaussian_kde
from scipy.optimize import fmin
from bisect import bisect_left, bisect_right
from scipy.stats import norm
import sys

class NoMaximaException(Exception):
    pass

def get_kde(data):
    kernel = gaussian_kde(data)
    kernel.silverman_factor()
    data_max = max(data)
    data_min = min(data)
    step = data_max / 300
    points = np.arange(data_min,data_max,step)
    return points,kernel.evaluate(points)

def half_gaussian(params,x):
    loc,scale,vol = params
    full_gauss = full_gaussian(params,x)
    half_point = bisect_left(x,loc)
    full_gauss[:half_point] = np.NaN
    return full_gauss

def full_gaussian(params,x):
    loc,scale,vol = params
    nd = norm(loc=loc,scale=scale)
    full_gauss = np.array(map(lambda q : nd.pdf(q)*vol,x))
    return full_gauss

def squared_difference(params,my_function,x,y):
    z = my_function(params,x)
    diff = [ d**2 for d in z-y ]
    return sum(diff)

def remove_zeros(data):
    cleaned_values = []
    for value in data:
        if value > 0.0:
            if not np.isfinite(value):
                print value
            cleaned_values.append(value)
    return cleaned_values

def log10_density(values):

    # Remove any genes with 0 FPKM
    cleaned_values = remove_zeros(values)

    # Get the log of the FPKM
    log_values = np.log10(cleaned_values)

    # Return the kernel density estimate
    return get_kde(log_values)

def get_derivative(y):
    dy = []
    lasty = y[0]
    for iy in y:
        dy.append( iy - lasty )
        lasty = iy
    return np.array(dy)

def get_maxima(x,y,stringency=2):
    dy = get_derivative(y)
    negative = [ iy < 0 for iy in dy ]
    positive = [ iy >= 0 for iy in dy ]
    points = []
    for i in range(len(x)):
        if i < stringency:
            continue
        if all(positive[i - stringency : i]) and all(negative[i : i + stringency]):
            points.append(i-1)
    if not points:
        raise NoMaximaException('Could not find any maxima with a stringency of %s')

    return points

def get_minima(x,y,stringency=2):
    dy = get_derivative(y)
    negative = [ iy < 0 for iy in dy ]
    positive = [ iy >= 0 for iy in dy ]
    points = []
    for i in range(len(x)):
        if i < stringency:
            continue
        if all(negative[i - stringency : i]) and all(positive[i : i + stringency]):
            points.append(i-1)
    return points

def get_last_maximum(x,y,stringency=2):

    maxima = get_maxima(x, y, stringency)

    return maxima[-1]

def recursive_last_maximum(x, y, stringency):
    if stringency == 0:
        raise NoMaximaException('Could not find any maxima at all!')
    try:
        last_max = get_last_maximum(x, y, stringency)
    except NoMaximaException:
        sys.stderr.write( 'couldnt get a maximum for stringency = {0}. Trying stringency = {1}\n'.format(stringency,stringency-1) )
        last_max = recursive_last_maximum(x, y, stringency-1)
    return last_max

def split_line(x, y, split_point):

    segment1_x = x[:split_point]
    segment1_y = y[:split_point]
    segment2_x = x[split_point:]
    segment2_y = y[split_point:]

    return (segment1_x, segment1_y) , (segment2_x, segment2_y)

def inverse_cumsum(y):
    cum_y = np.cumsum(y)
    max_y = max(cum_y)
    return max_y - cum_y

def remove_negatives( x ):
    for i in range(len(x)):
        if x[i] < 0:
            x[i] = 0
    return x

def sanitize_fdr(fdr):
    min_fdr = min(fdr)
    past_min = False
    for i in range(len(fdr)):
        if past_min:
            fdr[i] = 0
        if fdr[i] == min_fdr:
            past_min = True
    return fdr

def calculate_fdr(target_population, total):
    non_targets = remove_negatives( total - target_population )

    cumulative_non_targets = inverse_cumsum(non_targets)
    cum_targets = inverse_cumsum(target_population)

    fdr = cumulative_non_targets / cum_targets

    return sanitize_fdr(fdr)
