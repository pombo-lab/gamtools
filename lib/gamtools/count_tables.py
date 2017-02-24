"""
=======================
The count_tables module
=======================

The count_tables module contains non-cythonized (i.e. slow) functions for counts tables
(otherwise known as contingency tables). These functions are the equivalent
of those in the cosegregation_internal cython module, but they generally work for
n-dimensions and are orders of magnitude slower.

"""


import numpy as np


def get_transpositions(array):
    """
    Generator that iterates through all possible transpositions of
    an n-dimensional array.
    """

    axes = range(len(array.shape))
    for i in axes:
        yield tuple(axes[i:] + axes[:i])


def frequency_to_probability(counts_table):
    """
    Convert a contingency table expressed in frequencies to one
    expressed in probabilities.
    """

    total = counts_table.sum()
    probs_table = counts_table / float(total)

    return probs_table


def get_marginal_probabilities(probs_table):
    """
    Get the marginal probability of each event given a
    contingency table.
    """

    ind = []
    for transp in get_transpositions(probs_table):
        marginal_probs = [probs_table.transpose(transp)[1, ...].sum(),
                          probs_table.transpose(transp)[0, ...].sum()]
        ind.append(marginal_probs)
    return np.array(ind)


def either_locus_not_detected(probs):
    """
    Returns True if the probability of any event in a contingency
    table is 0.
    """

    return bool(probs.min())


def cosegregation(counts_table):
    """
    Return the co-segregation frequency of n loci given their
    contingency table.
    """

    probs_table = frequency_to_probability(counts_table)

    if either_locus_not_detected(probs_table):
        return np.NAN

    return probs_table.flat[-1]


def expected(counts_table):
    """
    Return the expected co-segregation probability of an arbitrary number
    of loci given their contingency table.
    """

    probs_table = frequency_to_probability(counts_table)
    marginal_probs = get_marginal_probabilities(probs_table)

    if either_locus_not_detected(marginal_probs):
        return np.NAN

    exp_freqs = marginal_probs.prod(axis=0)[0]
    return exp_freqs


def linkage(counts_table):
    """
    Return the linkage disequilibrium (D) for an arbitrary number of
    loci given their contingency table.
    """

    probs_table = frequency_to_probability(counts_table)
    marginal_probs = get_marginal_probabilities(probs_table)

    if either_locus_not_detected(marginal_probs):
        return np.NAN

    exp_freqs = marginal_probs.prod(axis=0)[0]
    observed = probs_table.flat[-1]
    if observed == 0:
        return np.NAN
    return observed - exp_freqs
