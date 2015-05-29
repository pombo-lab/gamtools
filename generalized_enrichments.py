import numpy as np
import itertools
import uuid
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Extract from a segmentation file only'
                                 'the positive windows for a particular sample')
parser.add_argument('-i','--interactions-file', required=True, 
                    help='A file containing all pairwise interactions between genomic windows')
parser.add_argument('-c','--classes-file', required=True, 
                    help='A file containing classes assigned to each genomic window')
parser.add_argument('-o','--output-filename', default='permutation_enrichment_results', 
                    help='First part of the output file name')
perm_type = parser.add_mutually_exclusive_group()
perm_type.add_argument('-p','--permutations', default=10, type=int, 
                    help='Number of times to randomly permute the input file')
perm_type.add_argument('-n','--no-permute', action='store_true',
                    help='Do not permute the input file at all')

def get_overlap(features_df1, features_df2, doublets_df):
    
    intermediate12 = pd.merge(features_df1, doublets_df, 
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_A'])
    overlap12 =       pd.merge(features_df2, intermediate12, 
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_B'])
    
    intermediate21 = pd.merge(features_df2, doublets_df, 
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_A'])
    overlap21 =       pd.merge(features_df1, intermediate21, 
                             left_on=['chrom', 'i'], right_on=['chrom', 'Pos_B'])
    
    both_overlap = pd.concat([overlap12[['chrom', 'Pos_A', 'Pos_B', 'Pi']],
                    overlap21[['chrom', 'Pos_A', 'Pos_B', 'Pi']]]).drop_duplicates()
    
    return both_overlap

def randomize_pos(pos, i, chrom_len):
    new_pos = pos + i
    if new_pos >= chrom_len:
        new_pos -= chrom_len
    return new_pos

def randomize_doublets(doublets_orig):
    
    doublets_df = doublets_orig.copy()
    
    i = np.random.randint(1, max_len)
    
    def randomize_chrom(chrom, i):
        
        chrom_len = chrom_lens[chrom]
        chrom_i = i % chrom_len

        doublets_df.ix[doublets_df.chrom == chrom,
                        ['Pos_A', 'Pos_B']] = np.array(
                            doublets_df.ix[doublets_df.chrom == chrom,
                                        ['Pos_A', 'Pos_B']]) + chrom_i
        
        doublets_df.ix[(doublets_df.chrom == chrom) & 
                       (doublets_df.Pos_A >= chrom_lens[chrom]),
                        ['Pos_A']] = np.array(doublets_df.ix[(doublets_df.chrom == chrom) & 
                                             (doublets_df.Pos_A >= chrom_lens[chrom]),
                                             ['Pos_A']]) - chrom_lens[chrom]
        
        doublets_df.ix[(doublets_df.chrom == chrom) & 
                       (doublets_df.Pos_B >= chrom_lens[chrom]),
                        ['Pos_B']] = np.array(doublets_df.ix[(doublets_df.chrom == chrom) & 
                                             (doublets_df.Pos_B >= chrom_lens[chrom]),
                                             ['Pos_B']]) - chrom_lens[chrom]
    
    for chrom in chroms:
        
        randomize_chrom(chrom, i)
        
    return doublets_df

def feature_pair_values(pairs_df, feat1, feat2):
    
    return get_overlap(window_classes[window_classes[feat1]],
                                   window_classes[window_classes[feat2]],
                                   pairs_df)

def get_feature_summary(pairs_df):
    
    results = []
    
    for feat1, feat2 in itertools.combinations_with_replacement(list(window_classes.columns[3:-1]), 2):
        
        feat_values = feature_pair_values(pairs_df, feat1, feat2)
        
        results.append((feat1, feat2, len(feat_values)))

    return results

def feature_summary(pairs_df, randomization_count):
    
    feature_pairs = list(itertools.combinations_with_replacement(list(window_classes.columns[3:-1]), 2))
    
    results = []
        
    for _r in range(randomization_count):

        rand_pairs = randomize_doublets(pairs_df)
        
        results.extend(get_feature_summary(rand_pairs))
        
    return results

def get_p_val(obs, rand):
    
    randoms = np.array(rand)
    
    p_high = 1.0 - (randoms > obs).mean()
    p_low = 1.0 - (randoms < obs).mean()

    return min((p_high, p_low))

if __name__ == '__main__':

    args = parser.parse_args()

    pairwise_interactions = pd.read_csv(args.interactions_file, delim_whitespace=True)

    window_classes = pd.read_csv(args.classes_file)

    chroms = ['chr{0}'.format(c) for c in range(1,20)]
    chrom_lens = window_classes.groupby('chrom').max().i.dropna()
    max_len = max(chrom_lens)

    if args.no_permute:
        results = get_feature_summary(pairwise_interactions)
    else:
        results = feature_summary(pairwise_interactions, args.permutations)

    results_df = pd.DataFrame.from_records(results, columns=['class1', 'class2', 'count'])

    if args.no_permute:
        results_df['permuted'] = 'no'
    else:
        results_df['permuted'] = 'yes'

    small_hash = str(uuid.uuid4())[:10]

    full_outpath = '{0}.{1}.csv'.format(args.output_filename, small_hash)

    results_df.to_csv(full_outpath, index=False)
