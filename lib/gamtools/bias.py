import argparse
import itertools

import pandas as pd
import numpy as np
from . import matrix
from mirnylib_numutils import removeDiagonals, observedOverExpected


def open_feature_file(feature_file_path):

    feature_bed = pd.read_csv(
        feature_file_path,
        sep='\t',
        names=['chrom', 'start', 'stop', 'feature']
    ).set_index(['chrom', 'start', 'stop'])

    return feature_bed


def calculate_bins(feature_bed):

    finite_vals = feature_bed.loc[
        np.isfinite(feature_bed.feature), 'feature'
    ]

    return ([-np.inf] + 
            [np.percentile(finite_vals, q) for q in range(10, 100, 10)] +
            [np.inf])


def discretize_feature(feature_bed, bins=10, percentiles=None):

    col_all_vals = np.array(feature_bed.feature).astype(np.float64)

    col_values = col_all_vals[np.isfinite(col_all_vals)]

    if not percentiles:
        percentiles = np.linspace(10,90,bins-1)

    col_bins = ([-np.inf] + 
                [np.percentile(col_values, q) for q in percentiles] +
                [np.inf])

    labels = ['bin_{}_to_{}'.format(low, high)
              for low, high
              in zip(col_bins[:-1], col_bins[1:])]

    return pd.cut(col_all_vals, bins=col_bins,
                  labels=labels)


def get_biases_heatmap(df, col, value_getter):

    col_labels = sorted([label for label in df.loc[~df[col].isnull(), col].unique()])
    bin_combinations = {get_key(l1,l2):[] for l1, l2 in itertools.combinations_with_replacement(col_labels, 2)}

    for chrom in ['chr{}'.format(c) for c in range(1,20)]:
        chrom_matrix = value_getter(chrom)

        for lab1, lab2 in bin_combinations.keys():

            bins1 = np.array(df.loc[df.chrom==chrom, col] == lab1)
            bins2 = np.array(df.loc[df.chrom==chrom, col] == lab2)

            bin_combinations[(lab1, lab2)].append(chrom_matrix[bins1,:][:,bins2])

    bin_comb_means = {labels:np.nanmean(np.concatenate([v.flatten() for v in values]))
                      for labels, values in bin_combinations.items()}

    n_cols = len(col_labels)
    means_matrix = np.zeros((n_cols, n_cols))
    for p,q in itertools.product(range(n_cols), range(n_cols)):
        means_matrix[p,q] = bin_comb_means[get_key(col_labels[p], col_labels[q])]

    return col_labels, means_matrix


def get_key(*vals):
    return tuple(sorted(vals))


def remove_zeros(a):
    s = np.nansum(a, axis=0) > 0
    b = a[:, s]
    c = b[s, :]
    return c, s


def replace_zeros(a, s, value=np.NAN):
    N = len(s)
    new_a = np.ones((N, N), dtype=a.dtype) * value
    tmp = np.ones((N, len(a)), dtype=a.dtype) * value
    tmp[s, :] = a
    new_a[:, s] = tmp
    return new_a


def get_obs_over_exp(A):

    removeDiagonals(A, 1)
    A, mask = remove_zeros(A)

    M = len(A.flat)
    toclip = 100 * min(0.999, (M - 10.) / M)

    A = observedOverExpected(A)
    A = np.clip(A, -1e10, np.percentile(A, toclip))

    return replace_zeros(A, mask)


def calculate_bias_matrix(feature_bed, matrices):

    col_labels = sorted([label for label in feature_bed.feature_bins.unique()])

    bin_combinations = {
        get_key(l1, l2):[] 
        for l1, l2 
        in itertools.combinations_with_replacement(col_labels, 2)
    }

    for matrix_path in matrices:
        (w1, w2), mat = matrix.read_file(matrix_path)

        normed_mat = get_obs_over_exp(mat)

        chrom_bins = feature_bed.loc[w1, 'feature_bins']

        for lab1, lab2 in bin_combinations.keys():

            bins1 = np.array(chrom_bins == lab1)
            bins2 = np.array(chrom_bins == lab2)

            bin_combinations[(lab1, lab2)].append(
                normed_mat[bins1,:][:,bins2])

    bin_comb_means = {labels:np.nanmean(np.concatenate(values, axis=None))
                      for labels, values in bin_combinations.items()}

    n_cols = len(col_labels)
    means_matrix = np.zeros((n_cols, n_cols))
    for p,q in itertools.product(range(n_cols), range(n_cols)):
        means_matrix[p,q] = bin_comb_means[get_key(col_labels[p], col_labels[q])]

    means_dataframe = pd.DataFrame(means_matrix)

    means_dataframe.index = col_labels
    means_dataframe.columns = col_labels

    return means_dataframe


def calculate_bias_from_args(args):

    feature_bed = open_feature_file(args.feature_path)

    feature_bed['feature_bins'] = discretize_feature(feature_bed)

    bias_matrix = calculate_bias_matrix(feature_bed, args.matrix_paths)

    bias_matrix.to_csv(args.output_path, sep='\t')
