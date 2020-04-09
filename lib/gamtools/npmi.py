import warnings

import pandas as pd
import numpy as np

def npmi_basic(A,B):
    freq_A = sum(A)/len(A)
    freq_B = sum(B)/len(B)
    freq_AB = np.logical_and(A, B).sum()/len(A)
    if not freq_AB:
        return np.NaN
    pmi = np.log2(freq_AB/(freq_A * freq_B))
    npmi = pmi/-np.log2(freq_AB)
    return(npmi)

def npmi_2d_slow(region1, region2):
    npmi_matrix = np.zeros((len(region1), len(region2))) * np.NaN
    for i1, window_1_arr in enumerate(region1):
        for i2, window_2_arr in enumerate(region2):
            npmi_matrix[i1, i2] = npmi_basic(window_1_arr, window_2_arr)
    return npmi_matrix

def npmi_2d_fast(region1, region2):
    M = region1.shape[1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pXY = region1.dot(region2.T) / M
        pXpY = (region1.sum(1).reshape(-1,1) / M) * (region2.sum(1) / M)
        pmi = np.log2(pXY / pXpY)
        npmi = pmi / -np.log2(pXY)
    return npmi

npmi_2d = npmi_2d_fast
