def get_segregation (coverage_per_bin, threshold):
    segregation_per_bin = coverage_per_bin > threshold
    segregation_per_bin = segregation_per_bin.astype(int)
    return segregation_per_bin

def get_threshold_for_sample (df, sample, percentile):
    import pandas as pd
    import numpy as np
    coverage_per_bin = df [sample]
    coverage_above0 = coverage_per_bin[coverage_per_bin > 75].tolist()
    coverage_per_bin_log10 = np.log10(coverage_above0)
    threshold = np.percentile (coverage_per_bin_log10, percentile)
    threshold_in_nucleotides = pow (10, threshold)
    return threshold_in_nucleotides

def get_percentile_based_segregation_table (coverage_df, threshold_set):
    import pandas as pd
    import numpy as np
    threshold_percentile = 100 - threshold_set
    df_segregation = coverage_df[coverage_df.columns[0:3]]
    list_of_samples = coverage_df.columns.tolist()
    list_of_samples = list_of_samples[3:]
    list_of_thresholds = []
    for sample in list_of_samples:
        coverage_per_bin = coverage_df [sample]
        coverage_above0 = coverage_per_bin[coverage_per_bin > 75].tolist()
        coverage_per_bin_log10 = np.log10(coverage_above0)
        threshold = np.percentile (coverage_per_bin_log10, threshold_percentile)
        threshold_in_nucleotides = pow (10, threshold)
        sample_segregation = get_segregation (coverage_per_bin, threshold_in_nucleotides)
        df_segregation = pd.concat([df_segregation, sample_segregation], axis=1, join_axes=[df_segregation.index])
    return df_segregation

def get_percentile_threshold_for_dataset (reads_df, coverage_df):
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    from random import shuffle 
    plt.switch_backend('agg')
    list_of_samples = coverage_df.columns.tolist()
    if len (list_of_samples) < 303:
        list_of_samples = list_of_samples [3:] 
    if len (list_of_samples) > 303:
        list_of_samples = list_of_samples [3:303] 
        shuffle (list_of_samples)
    percentages_over_bins_mean, percentages_over_bins_std = [], []
    #bins_to_calculate = [6, 9] # Correspond to 2 and 3 reads per bin
    bins_to_calculate = [6] # Correspond to 2 reads per bin
    print ('Starting to calculate threshold percentile...')
    for bar in bins_to_calculate:
        percentages_over_samples_mean, percentage_over_samples_stdev = [], []
        percentages_over_samples = []
        for sample in list_of_samples:
            sys.stdout.write(".")
            sys.stdout.flush()
            coverage_per_sample = coverage_df [sample].tolist()
            reads_per_sample = reads_df [sample].tolist()
            percentages_over_percentile = []
            for percentile in reversed (range (1,100)):
                threshold_for_sample = get_threshold_for_sample (coverage_df, sample, percentile)
                reads_in_all_windows, reads_in_positive_windows = [], []
                for reads, coverage in zip(reads_per_sample, coverage_per_sample):
                    if coverage > 0: # Take away windows with no reads
                        reads_in_all_windows.append(reads)
                    if coverage >= threshold_for_sample:# positive windows
                        reads_in_positive_windows.append(reads)
                all_win = np.log10(reads_in_all_windows)
                pos_win = np.log10(reads_in_positive_windows)
                bins_all = np.histogram(np.hstack((all_win,pos_win)), bins=80, range = (0, 4), density=False)[1]
                counts_all_list = plt.hist(all_win, bins_all)[0]
                counts_pos_list = plt.hist(pos_win, bins_all)[0]
                if counts_all_list [bar] == 0:
                    percent_bar = 0
                else:
                    percent_bar = ((counts_pos_list [bar] / counts_all_list [bar]) * 100)
                plt.close('all')
                percentages_over_percentile.append (percent_bar)
            percentages_over_samples.append (percentages_over_percentile)
        percentages_over_samples_mean = np.mean (percentages_over_samples, axis=0)
        percentage_over_samples_stdev = np.std (percentages_over_samples, axis=0)
        percentages_over_bins_mean.append (percentages_over_samples_mean) # for more than one bin
        percentages_over_bins_std.append (percentage_over_samples_stdev)
    for mean, k in zip(percentages_over_bins_mean [0], range(0, len(percentages_over_bins_mean [0]))):
        if mean > 0.1:
            threshold_set = k
            break
    return threshold_set, percentages_over_bins_mean, percentages_over_bins_std

def GAM_segregation_tables (reads_df, coverage_df, name_of_segregation_table):
    threshold_set, percentages_over_bins_mean, percentages_over_bins_std = get_percentile_threshold_for_dataset (reads_df, coverage_df)
    segregation_table = get_percentile_based_segregation_table (coverage_df, threshold_set)
    name_of_segregation_table = name_of_segregation_table + '.tp' + str (threshold_set) + '.table'
    segregation_table.to_csv (name_of_segregation_table, sep='\t', index=False, header=True)
    print ('Threshold percentile for the dataset is ' + str(threshold_set))
    print ('Segregation table was saved as ' + name_of_segregation_table)


