"""
=====================
The enrichment module
=====================

Identifying patterns in chromatin interaction datasets
======================================================

The enrichment module contains functions for identifying sets
of genomic features which are over- or under-represented in
chromatin interactions. In essence, given a list of interactions
of interest between genomic windows, and a classification of those
same windows (for example, windows could be classified according
to whether or not they contain CTCF binding sites),
:func:`do_enrichment` function can be used to randomly permute
the input interactions and determine whether interactions between
different classes of window are more or less likely than expected
by chance.

"""

import uuid
import itertools

import numpy as np
import pandas as pd


def get_overlap(features_df1, features_df2, doublets_df):
    """Subset a table of interactions given two tables of possible interacting elements.

    Takes a table of pairwise interactions (doublets_df) and returns
    only those interactions which involve a window from
    features_df1 interacting with a window from features_df2.

    :param features_df1: Table of windows allowed on one side of the interaction \
            (can be either left or right).
    :type features_df1: :class:`pandas.DataFrame`
    :param features_df2: Table of windows allowed on the other side of the interaction.
    :type features_df2: :class:`pandas.DataFrame`
    :param doublets_df: Table of pairwise interactions.
    :type doublets_df: :class:`pandas.DataFrame`
    :returns: A subset of the original doublets_df.

    >>> pairwise_interactions = pd.DataFrame([('chr1', 10, 20, 0.75),
    ...                                        ('chr2', 10, 20, 0.5),
    ...                                        ('chr1', 10, 30, 0.4)],
    ...                                      columns=['chrom', 'Pos_A', 'Pos_B',
    ...                                               'interaction'])
    >>> pairwise_interactions
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75
    1  chr2     10     20         0.50
    2  chr1     10     30         0.40
    >>> enhancers = pd.DataFrame([('chr1', 10),
    ...                           ('chr2', 20)],
    ...                          columns=['chrom', 'i'])
    >>> enhancers
      chrom   i
    0  chr1  10
    1  chr2  20
    >>> genes = pd.DataFrame([('chr1', 20),
    ...                       ('chr2', 30)],
    ...                      columns=['chrom', 'i'])
    >>> genes
      chrom   i
    0  chr1  20
    1  chr2  30
    >>> enrichment.get_overlap(enhancers, genes, pairwise_interactions)
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75

    """

    intermediate12 = pd.merge(
        features_df1, doublets_df, left_on=[
            'chrom', 'i'], right_on=[
                'chrom', 'Pos_A'])
    overlap12 = pd.merge(features_df2, intermediate12,
                         left_on=['chrom', 'i'], right_on=['chrom', 'Pos_B'])

    intermediate21 = pd.merge(
        features_df2, doublets_df, left_on=[
            'chrom', 'i'], right_on=[
                'chrom', 'Pos_A'])
    overlap21 = pd.merge(features_df1, intermediate21,
                         left_on=['chrom', 'i'], right_on=['chrom', 'Pos_B'])

    both_overlap = pd.concat([overlap12[['chrom', 'Pos_A', 'Pos_B', 'interaction']],
                              overlap21[['chrom', 'Pos_A', 'Pos_B', 'interaction']]
                             ]).drop_duplicates()

    return both_overlap


def randomize_doublets(doublets_orig, chrom_lengths):
    """Randomize a table of pairwise interactions.

    Randomize a table of pairwise interactions by shifting all interactions
    on a given chromosome by the same amount. This preserves both
    the span of each individual interaction and the distance distribution
    between different interactions. Shifts are done by circular
    permutation, so if either end of an interaction would be shifted
    off the end of the chromosome, it will be replaced at the
    beginning. This does mean that a small number of interactions
    might be split in half by the randomization (i.e. the end is
    moved to the beginning of the chromosome).

    :param doublets_orig: Table of pairwise interactions.
    :type doublets_orig: :class:`pandas.DataFrame`
    :param dict chrom_lengths: Dictionary giving the number of windows \
            in each chromosome.
    :returns: Randomized version of doublets_orig.

    >>> pairwise_interactions = pd.DataFrame([('chr1', 10, 20, 0.75),
    ...                                        ('chr2', 10, 20, 0.5),
    ...                                        ('chr1', 10, 30, 0.4)],
    ...                                      columns=['chrom', 'Pos_A', 'Pos_B',
    ...                                               'interaction'])
    >>> pairwise_interactions
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75
    1  chr2     10     20         0.50
    2  chr1     10     30         0.40
    >>> chromosome_lengths = {'chr1':50, 'chr2':35}
    >>> enrichment.randomize_doublets(pairwise_interactions, chromosome_lengths)
      chrom  Pos_A  Pos_B  interaction
    0  chr1     38     48         0.75
    1  chr2      3     13         0.50
    2  chr1     38      8         0.40

    """

    # TODO: Handle split interactions where one half is permuted past the
    # chromosome length

    doublets_df = doublets_orig.copy()

    shift = np.random.randint(1, max(chrom_lengths.values()))

    def randomize_chrom(chrom, shift):
        """
        Internal function that shifts all interactions on a given
        chromosome `chrom` by an amount `shift`.
        """

        chrom_len = chrom_lengths[chrom]

        # If a shift is bigger than the length of the chromosome,
        # the end result is the same as shifting by the remainder
        # when shift is divided by the chromosome length
        chrom_shift = shift % chrom_len

        # Shift both positions by chrom_shift
        doublets_df.ix[doublets_df.chrom == chrom,
                       ['Pos_A', 'Pos_B']] = np.array(
                           doublets_df.ix[doublets_df.chrom == chrom,
                                          ['Pos_A', 'Pos_B']]) + chrom_shift

        # Check left doesn't extend past chromosome end
        doublets_df.ix[
            (doublets_df.chrom == chrom) & (
                doublets_df.Pos_A >= chrom_lengths[chrom]),
            ['Pos_A']] = np.array(
                doublets_df.ix[
                    (doublets_df.chrom == chrom) & (
                        doublets_df.Pos_A >= chrom_lengths[chrom]),
                    ['Pos_A']]) - chrom_lengths[chrom]

        # Check right doesn't extend past chromosome end
        doublets_df.ix[
            (doublets_df.chrom == chrom) & (
                doublets_df.Pos_B >= chrom_lengths[chrom]),
            ['Pos_B']] = np.array(
                doublets_df.ix[
                    (doublets_df.chrom == chrom) & (
                        doublets_df.Pos_B >= chrom_lengths[chrom]),
                    ['Pos_B']]) - chrom_lengths[chrom]

    for chrom in chrom_lengths.keys():

        randomize_chrom(chrom, shift)

    return doublets_df


def feature_pair_values(pairwise_interactions, window_classes, feat1, feat2):
    """Subset a table of interactions, to obtain only those where windows
    from the first classification interact with those from the second classification.

    Takes a table of pairwise interactions (doublets_df), a table of window
    classifications and two specific window classifications. Returns only those
    interactions which involve a window classified as feat1 interacting with a
    window classified as feat2.

    :param pairwise_interactions: Table of pairwise interactions.
    :type pairwise_interactions: :class:`pandas.DataFrame`
    :param str feat1: Classification of windows allowed on one side of the interaction \
            (can be either left or right).
    :param str feat2: Classification of windows allowed on the other side of the interaction.
    :returns: A subset of the original doublets_df.

    >>> pairwise_interactions = pd.DataFrame([('chr1', 10, 20, 0.75),
    ...                                        ('chr2', 10, 20, 0.5),
    ...                                        ('chr1', 10, 30, 0.4)],
    ...                                      columns=['chrom', 'Pos_A', 'Pos_B',
    ...                                               'interaction'])
    >>> pairwise_interactions
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75
    1  chr2     10     20         0.50
    2  chr1     10     30         0.40
    >>> window_classification = pd.DataFrame([('chr1', 10, True, False),
    ...                                       ('chr1', 20, False, True),
    ...                                       ('chr2', 20, True, False),
    ...                                       ('chr2', 30, False, True)],
    ...                                      columns=['chrom', 'i', 'Enhancer', 'Gene'])
    >>> window_classification
      chrom   i Enhancer   Gene
    0  chr1  10     True  False
    1  chr1  20    False   True
    2  chr2  20     True  False
    3  chr2  30    False   True
    >>> enrichment.feature_pair_values(pairwise_interactions, window_classification,
    ...                                'Enhancer', 'Gene')
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75

    """

    return get_overlap(window_classes[window_classes[feat1]],
                       window_classes[window_classes[feat2]],
                       pairwise_interactions)


def get_feature_summary(pairwise_interactions, window_classes):
    """
    Count the number of pairwise interactions between different window classes.

    :param pairwise_interactions: Table of pairwise interactions.
    :type pairwise_interactions: :class:`pandas.DataFrame`
    :param window_classes: Table of windows and their classifications
    :type window_classes: :class:`pandas.DataFrame`
    :returns: A list of tuples of the form (first class, second class, count).

    >>> pairwise_interactions = pd.DataFrame([('chr1', 10, 20, 0.75),
    ...                                        ('chr2', 10, 20, 0.5),
    ...                                        ('chr1', 10, 30, 0.4)],
    ...                                      columns=['chrom', 'Pos_A', 'Pos_B',
    ...                                               'interaction'])
    >>> pairwise_interactions
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75
    1  chr2     10     20         0.50
    2  chr1     10     30         0.40
    >>> window_classification = pd.DataFrame([('chr1', 10, True, False),
    ...                                       ('chr1', 20, False, True),
    ...                                       ('chr2', 20, True, False),
    ...                                       ('chr2', 30, False, True)],
    ...                                      columns=['chrom', 'i', 'Enhancer', 'Gene'])
    >>> window_classification
      chrom   i Enhancer   Gene
    0  chr1  10     True  False
    1  chr1  20    False   True
    2  chr2  20     True  False
    3  chr2  30    False   True
    >>> enrichment.get_feature_summary(pairwise_interactions, window_classification)
    [('Enhancer', 'Enhancer', 0), ('Enhancer', 'Gene', 1), ('Gene', 'Gene', 0)]
    """

    results = []
    feature_classes = [col for col in window_classes.columns
                       if not col in ['chrom', 'start', 'stop', 'i']]

    for feat1, feat2 in itertools.combinations_with_replacement(
            feature_classes, 2):

        feat_values = feature_pair_values(
            pairwise_interactions, window_classes, feat1, feat2)

        results.append((feat1, feat2, len(feat_values)))

    return results


def randomized_summary(
        pairwise_interactions,
        window_classes,
        chrom_lengths,
        randomization_count):
    """
    Randomly shuffle a table of pairwise interactions and count the number of
    randomized interactions involving each feature class.

    :param pairwise_interactions: Table of pairwise interactions.
    :type pairwise_interactions: :class:`pandas.DataFrame`
    :param window_classes: Table of windows and their classifications
    :type window_classes: :class:`pandas.DataFrame`
    :param dict chrom_lengths: Dictionary giving the number of windows \
            in each chromosome.
    :param int randomization_count: Number of times to randomize.
    :returns: A list of tuples of the form (first class, second class, count).

    >>> pairwise_interactions = pd.DataFrame([('chr1', 10, 20, 0.75),
    ...                                        ('chr2', 10, 20, 0.5),
    ...                                        ('chr1', 10, 30, 0.4)],
    ...                                      columns=['chrom', 'Pos_A', 'Pos_B',
    ...                                               'interaction'])
    >>> pairwise_interactions
      chrom  Pos_A  Pos_B  interaction
    0  chr1     10     20         0.75
    1  chr2     10     20         0.50
    2  chr1     10     30         0.40
    >>> window_classification = pd.DataFrame([('chr1', 10, True, False),
    ...                                       ('chr1', 20, False, True),
    ...                                       ('chr2', 20, True, False),
    ...                                       ('chr2', 30, False, True)],
    ...                                      columns=['chrom', 'i', 'Enhancer', 'Gene'])
    >>> window_classification
      chrom   i Enhancer   Gene
    0  chr1  10     True  False
    1  chr1  20    False   True
    2  chr2  20     True  False
    3  chr2  30    False   True
    >>> chromosome_lengths = {'chr1':50, 'chr2':35}
    >>> enrichment.randomized_summary(pairwise_interactions, window_classification,
    ...                               chromosome_lengths, 2)
    [('Enhancer', 'Enhancer', 0),
     ('Enhancer', 'Gene', 0),
     ('Gene', 'Gene', 0),
     ('Enhancer', 'Enhancer', 0),
     ('Enhancer', 'Gene', 0),
     ('Gene', 'Gene', 0)]
    """

    results = []

    for _ in range(randomization_count):

        rand_pairs = randomize_doublets(pairwise_interactions, chrom_lengths)

        results.extend(get_feature_summary(rand_pairs, window_classes))

    return results


def get_p_val(obs, rand):
    """Calculate a two-tailed p-value from an observed value and a
    list of permuted (bootstrapped) values.

    :param float obs: Observed value.
    :param list rand: List of values obtained by permutation.
    :returns: The probability of observing a value as large (or as \
            small) as obs given the background distribution rand.
    """

    randoms = np.array(rand)

    p_high = 1.0 - (randoms > obs).mean()
    p_low = 1.0 - (randoms < obs).mean()

    return min((p_high, p_low))


def get_full_output_path(output_prefix, num_permutations):
    """Generate the path for an output csv file.

    The path contains a random string so that permutations can
    be easily parallelized without overwriting each other.

    :param str output_prefix: Beginning part of the csv file path.
    :param int num_permutations: Number of times the data was permuted.
    :returns: Output path in the format <prefix>_<num_permutations>.<random_string>.csv

    >>> enrichment.get_full_output_path('gam_gene_enrichments', 50)
    'gam_gene_enrichments_50_permutations.8d73b9dd-2.csv'
    >>> enrichment.get_full_output_path('gam_gene_enrichments', 0)
    'gam_gene_enrichments_unpermuted.d059ad6c-2.csv'
    """

    if num_permutations == 0:
        permute_string = 'unpermuted'
    else:
        permute_string = '{}_permutations'.format(num_permutations)

    small_hash = str(uuid.uuid4())[:10]

    # TODO: Check that this file doesn't exist
    full_outpath = '{0}_{1}.{2}.csv'.format(
        output_prefix, permute_string, small_hash)

    return full_outpath


def do_enrichment(
        pairwise_interactions,
        window_classes,
        num_permutations,
        output_prefix,
        chroms=None):
    """Calculate the number of pairwise pairwise_interactions connecting each class of windows
    (either in real or permuted data) and save the results to a file

    :param pairwise_interactions: Table of pairwise interactions.
    :type pairwise_interactions: :class:`pandas.DataFrame`
    :param window_classes: Table of windows and their classifications
    :type window_classes: :class:`pandas.DataFrame`
    :param int num_permutations: Number of times to randomize the data.
    :param str output_prefix: Beginning part of the csv file path.
    :param list chroms: List of chromosomes to consider (default is \
            the mouse autosomes, i.e. chr1 to chr19).
    """

    if chroms is None:
        chroms = ['chr{0}'.format(c) for c in range(1, 20)]

    chrom_lengths = {
        chrom: len(
            window_classes[
                window_classes.chrom == chrom]) for chrom in chroms}

    if num_permutations == 0:
        results = get_feature_summary(pairwise_interactions, window_classes)
    else:
        results = randomized_summary(
            pairwise_interactions,
            window_classes,
            chrom_lengths,
            num_permutations)

    results_df = pd.DataFrame.from_records(
        results, columns=['class1', 'class2', 'count'])

    if num_permutations == 0:
        results_df['permuted'] = 'no'
    else:
        results_df['permuted'] = 'yes'

    full_outpath = get_full_output_path(output_prefix, num_permutations)

    results_df.to_csv(full_outpath, index=False)


def enrichment_from_args(args):
    """
    Wrapper function to call do_enrichment from argparse.
    """

    interactions = pd.read_csv(args.interactions_file, delim_whitespace=True)
    if 'interaction' not in interactions.columns:
        interactions = interactions.rename(columns={'Pi': 'interaction'})

    window_classes = pd.read_csv(args.classes_file)
    do_enrichment(interactions, window_classes,
                  args.num_permutations, args.output_prefix)
