import numpy as np

########################################################
#
# Non-cythonized (i.e. slow) functions for counts tables
#
########################################################

def get_transpositions(array):
    axes = range(len(array.shape))
    for i in axes:
        yield tuple(axes[i:] + axes[:i])


def frequency_to_probability(counts_table):

    total = counts_table.sum()
    probs_table = counts_table / float(total)

    return probs_table


def get_marginal_probabilities(probs_table):

    ind = []
    for t in get_transpositions(probs_table):
        marginal_probs = [probs_table.transpose(t)[1,...].sum(),
                 probs_table.transpose(t)[0,...].sum()]
        ind.append(marginal_probs)
    return np.array(ind)


def either_locus_not_detected(probs):

    if probs.min() == 0.0:
        return True
    else:
        return False


def cosegregation(counts_table):

    probs_table = frequency_to_probability(counts_table)

    if either_locus_not_detected(probs_table):
        return np.NAN

    return probs_table.flat[-1]

def expected(counts_table):
    probs_table = frequency_to_probability(counts_table)
    marginal_probs = get_marginal_probabilities(probs_table)

    if either_locus_not_detected(marginal_probs):
        return np.NAN

    expected = marginal_probs.prod(axis=0)[0]
    return expected


def D(counts_table):
    probs_table = frequency_to_probability(counts_table)
    marginal_probs = get_marginal_probabilities(probs_table)

    if either_locus_not_detected(marginal_probs):
        return np.NAN

    expected = marginal_probs.prod(axis=0)[0]
    observed = probs_table.flat[-1]
    if observed == 0:
        return np.NAN
    return observed - expected
