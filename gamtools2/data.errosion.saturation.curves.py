"""
Functions perform data errosion test on a segregation table and calcuate how many pairs 
of loci are seen at-least once across the dataset with different number of samples.
Performing the test at different resolutions suggests an optimal resolution for the current size of the dataset.
"""

def get_number_of_pairs_never_seen_together_for_region (segregation_table, region):
    from gamtools import cosegregation
    coseg = cosegregation.get_cosesgregation (segregation_table, region)
    coseg_df = pd.DataFrame (coseg)
    pairs_not_seen_together = (coseg_df == 0).values.sum()
    pairs_seen_at_least_once = (coseg_df > 0).values.sum()
    total_pairs = (coseg_df >= 0).values.sum()
    percent_never_seen = (pairs_not_seen_together / total_pairs)*100
    percent_never_seen = (round (percent_never_seen, 3))
    percent_seen_at_least_once = (pairs_seen_at_least_once / total_pairs)*100
    percent_seen_at_least_once = (round (percent_seen_at_least_once, 3))
    return percent_never_seen, percent_seen_at_least_once

def percent_of_pairs_never_seen_together_data_errosion (segregation_table, list_of_samples, start_samples, step):
    import pandas as pd
    from random import shuffle
    list_percent_seen_at_least_once = []
    number_of_samples = len(segregation_table.columns)
    myrange = range(start_samples, number_of_samples, step)
    for samples_in_table in myrange:
        subsampling = list_of_samples [0:samples_in_table]
        sub_table = segregation_table.loc[:,subsampling]
        percent_never_seen, percent_seen_at_least_once = get_number_of_pairs_never_seen_together_for_region (sub_table, 'chr19')
        list_percent_seen_at_least_once.append (percent_seen_at_least_once)
    return list_percent_seen_at_least_once

def zero_slope (data, chunksize=5, max_slope = 0.1, step=5):
    midindex = int(chunksize / 2)
    for index in range(len(data) - chunksize):
        chunk = data[index:index + chunksize]
        slope_change = (abs(chunk[-1] - chunk[0]))
        if 0 < slope_change < max_slope:
            saturation_point = str(index*step)
            saturation_point = saturation_point + ' samples'
            break
        if index == (len(data) - chunksize -1):
            saturation_point = 'never'
    return saturation_point





