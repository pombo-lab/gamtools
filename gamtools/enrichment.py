import uuid
import numpy as np
import pandas as pd
import itertools

def get_overlap(features_df1, features_df2, doublets_df):

    intermediate12 = pd.merge(features_df1, doublets_df,
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_A'])
    overlap12 =       pd.merge(features_df2, intermediate12,
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_B'])

    intermediate21 = pd.merge(features_df2, doublets_df,
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_A'])
    overlap21 =       pd.merge(features_df1, intermediate21,
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_B'])

    both_overlap = pd.concat([overlap12[['chrom', 'Pos_A', 'Pos_B', 'interaction']],
                    overlap21[['chrom', 'Pos_A', 'Pos_B', 'interaction']]]).drop_duplicates()

    return both_overlap

def randomize_pos(pos, i, chrom_len):
    new_pos = pos + i
    if new_pos >= chrom_len:
        new_pos -= chrom_len
    return new_pos

def randomize_doublets(doublets_orig, chrom_lengths):

    doublets_df = doublets_orig.copy()

    i = np.random.randint(1, max(chrom_lengths.values()))

    def randomize_chrom(chrom, i):

        chrom_len = chrom_lengths[chrom]
        chrom_i = i % chrom_len

        doublets_df.ix[doublets_df.chrom == chrom,
                        ['Pos_A', 'Pos_B']] = np.array(
                            doublets_df.ix[doublets_df.chrom == chrom,
                                        ['Pos_A', 'Pos_B']]) + chrom_i

        doublets_df.ix[(doublets_df.chrom == chrom) &
                       (doublets_df.Pos_A >= chrom_lengths[chrom]),
                        ['Pos_A']] = np.array(doublets_df.ix[(doublets_df.chrom == chrom) &
                                             (doublets_df.Pos_A >= chrom_lengths[chrom]),
                                             ['Pos_A']]) - chrom_lengths[chrom]

        doublets_df.ix[(doublets_df.chrom == chrom) &
                       (doublets_df.Pos_B >= chrom_lengths[chrom]),
                        ['Pos_B']] = np.array(doublets_df.ix[(doublets_df.chrom == chrom) &
                                             (doublets_df.Pos_B >= chrom_lengths[chrom]),
                                             ['Pos_B']]) - chrom_lengths[chrom]

    for chrom in chrom_lengths.keys():

        randomize_chrom(chrom, i)

    return doublets_df

def feature_pair_values(pairs_df, window_classes, feat1, feat2):

    return get_overlap(window_classes[window_classes[feat1]],
                       window_classes[window_classes[feat2]],
                                   pairs_df)

def get_feature_summary(pairs_df, window_classes):

    results = []

    for feat1, feat2 in itertools.combinations_with_replacement(list(window_classes.columns[3:-1]), 2):

        feat_values = feature_pair_values(pairs_df, window_classes, feat1, feat2)

        results.append((feat1, feat2, len(feat_values)))

    return results

def randomized_summary(pairs_df, window_classes, chrom_lengths, randomization_count):

    results = []

    for _r in range(randomization_count):

        rand_pairs = randomize_doublets(pairs_df, chrom_lengths)

        results.extend(get_feature_summary(rand_pairs, window_classes))

    return results

def get_p_val(obs, rand):

    randoms = np.array(rand)

    p_high = 1.0 - (randoms > obs).mean()
    p_low = 1.0 - (randoms < obs).mean()

    return min((p_high, p_low))

def get_full_output_path(output_prefix, num_permutations):

    if num_permutations == 0:
        permute_string = 'unpermuted'
    else:
        permute_string = '{}_permutations'.format(num_permutations)

    small_hash = str(uuid.uuid4())[:10]

    full_outpath = '{0}_{1}.{2}.csv'.format(output_prefix, permute_string, small_hash)

    return full_outpath

def do_enrichment(interactions, window_classes, num_permutations, output_prefix, chroms=None):

    if chroms is None:
        chroms = ['chr{0}'.format(c) for c in range(1,20)]

    chrom_lengths = {chrom:len(window_classes[window_classes.chrom == chrom]) for chrom in chroms}

    if num_permutations == 0:
        results = get_feature_summary(interactions, window_classes)
    else:
        results = randomized_summary(interactions, window_classes, chrom_lengths, num_permutations)

    results_df = pd.DataFrame.from_records(results, columns=['class1', 'class2', 'count'])

    if num_permutations == 0:
        results_df['permuted'] = 'no'
    else:
        results_df['permuted'] = 'yes'

    full_outpath = get_full_output_path(output_prefix, num_permutations)

    results_df.to_csv(full_outpath, index=False)


def enrichment_from_args(args):

    interactions = pd.read_csv(args.interactions_file, delim_whitespace=True)
    if 'interaction' not in interactions.columns:
        interactions = interactions.rename(columns={'Pi':'interaction'})

    window_classes = pd.read_csv(args.classes_file)
    do_enrichment(interactions, window_classes,
                  args.num_permutations, args.output_prefix)

