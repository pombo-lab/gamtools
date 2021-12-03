"""
======================
The permutation module
======================

The permutation module contains functions for randomly permuting GAM
:ref:`segregation tables <segregation_table>`, which can be a useful
way of generating random backgrounds for comparison with real datasets.

"""

import numpy as np
import pandas as pd

from . import segregation

def permute_by_offset(sample_segregation, offset):
    """Circularly permute a single column of a segregation table.

    This function takes one single column from a
    :ref:`segregation_table` (i.e. one sample, or one NP) and circularly
    permutes it by "offset" bins.

    :param sample_segregation: Input column of a segregation table to permute.
    :param int offset: Number of bins to permute by.
    :returns: Returns a newly randomized :ref:`segregation_table`
    """

    offset = offset % len(sample_segregation)
    # Moving each value in an array of length L right by x bins is the same
    # as splitting the data at bin (L - x) and swapping the two halves
    corrected_offset = len(sample_segregation) - offset
    new_start = sample_segregation.iloc[corrected_offset:]
    new_end = sample_segregation.iloc[:corrected_offset]
    new_col = pd.concat([new_start, new_end]).values

    return new_col


def permute_by_chromosome(sample_segregation, offset):
    """Separately permute each chromosome from a single column of a segregation table.

    This function takes one single column from a
    :ref:`segregation_table` (i.e. one sample, or one NP) and circularly
    permutes each chromosome separately by "offset" bins.

    :param sample_segregation: Input column of a segregation table to permute.
    :param int offset: Number of bins to permute by.
    :returns: Returns a newly randomized :ref:`segregation_table`
    """

    permuted_chromosomes = []
    chrom_index = sample_segregation.index.get_level_values(0)

    for chrom in chrom_index.unique():
        original_chromosome = sample_segregation[chrom_index == chrom]
        permuted_chromosome = permute_by_offset(original_chromosome, offset)
        permuted_chromosomes.append(permuted_chromosome)

    permuted_segregation = np.concatenate(permuted_chromosomes)

    return permuted_segregation


def permute_segregation(input_segregation):
    """Circularly permute each column of a segregation table.

    This function takes a table of GAM segregation data (see
    :ref:`segregation_table`) and circularly permutes each column by a random
    amount. For example, if the column contains [ 0, 0, 0, 1, 1, 1 ] then a
    circular permutation by one unit would give: [ 1, 0, 0, 0, 1, 1 ]. This can
    be useful because it preserves the scaling features of a GAM matrix (i.e.
    windows which lie next to each other are still more likely to co-segregate
    in the same tube) but it randomizes long-range interactions.

    Another feature of "real" GAM data is that it contains unmappable regions
    with no signal. In order to preserve this feature, genomic regions which
    are never detected in the input :ref:`segregation_table` are  not subjected
    to permutation.

    :param input_segregation: Input segregation table to permute.
    :type input_segregation: :ref:`segregation_table`
    :returns: Returns a newly randomized :ref:`segregation_table`
    """

    # Only permute positions that were mapped at least once
    # This means we preserve the locations of centromeres etc,
    # but requires that the input data has a large number of
    # columns
    mappable = input_segregation.sum(axis=1).astype(bool)
    no_windows, no_samples = input_segregation[mappable].shape

    # Make a copy of the original data
    permutation = input_segregation.copy()

    # Loop over columns
    for i in range(no_samples):

        # Choose a random position to break the segregation in two,
        # Swap the two chunks around and write them to the copied df
        offset = np.random.randint(no_windows)

        new_col = permute_by_chromosome(input_segregation.iloc[mappable.values, i], offset)

        permutation.iloc[mappable.values, i] = new_col

    return permutation


def permute_segregation_autosomal(input_segregation, autosomes=None):
    """Circularly permute each autosomal chromosome in a segregation table

    This function takes a table of GAM segregation data (see
    :ref:`segregation_table`) and circularly permutes each autosomal chromosome
    by a random amount. For each GAM sample, The last n values covering each
    chromosome are moved from the end of the chromosome to the beginning, where
    n is a random integer between 0 and the length of the chromosome.

    Separately permuting each chromosome preserves their individual detection
    frequencies, and avoids permuting genomic regions with very different
    detection frequencies (e.g. mitochondrial DNA and sex chromosomes) into
    one another.

    :param input_segregation: Input segregation table to permute.
    :type input_segregation: :ref:`segregation_table`
    :returns: Returns a newly randomized :ref:`segregation_table`
    """

    # Sex chromosomes, unjoined contigs and mtDNA behave weirdly,
    # so don't permute them into the autosomal regions.

    autosomes = ['chr{0}'.format(c) for c in range(1, 20)]

    is_autosomal = input_segregation.index.get_level_values(0).isin(autosomes)

    segregation_to_permute = input_segregation.loc[is_autosomal]

    permuted_segregation = permute_segregation(segregation_to_permute)

    return permuted_segregation


def permute_segregation_from_args(args):
    """Extract parameters from an argparse namespace object and pass them to
    permute_segregation_autosomal.
    """

    input_segregation = segregation.open_segregation(args.segregation_file)

    permuted_segregation = permute_segregation_autosomal(input_segregation)

    permuted_segregation.to_csv(
        args.output_file,
        sep='\t',
        header=True,
        index=True)
