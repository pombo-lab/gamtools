'''''''''''
Script by Tom 
19/10/2018

This script is to take segregation tables from GAM data and output metrics to 
see whether there is potential contamination

similarity between wells is measured in Jaccard distance and clustered using the default 
seaborn clustermap settings

**** WRAPPER FUNCTION IS AT THE BOTTOM OF THIS SCRIPT ****

The only neccessary input is the name and location of a segregation table containing 
all and only samples for that plate

If your file contains more than one plate, use the plate_identifier 
option in the function and supply a plate identity as a string. The plate identity is a
short string that is present in the names of all samples that belong to one plate

**** options for wrapper function: ****

input_table = location and file of segregation table
plate_identifier = None or string identifier for plate
similarity_threshold = default is 0.1 , i.e return all wells which have over 90% similarity
dpi_value = default is 500, qualtiy for plot produces
fig_output_destination ='outputfig.png', name of file for the sns figure
bad_well_file_destination = 'list_of_bad_wells', name of file for list of wells that failed the threshold as txt file
output_raw_df = 'raw_values.txt', name of file for raw values of the clustered df 

example:
contamination_bias_wrapper('your_segregation_table',plate_identifier = 'plate_identity_1')
'''''''''''

import numpy as np
import pandas as pd
import scipy.spatial
jaccard = scipy.spatial.distance.jaccard
import matplotlib.pyplot as plt
import seaborn as sns

# function for subsetting out columns of dataframe that contain a certain string
def get_plate(data_segregation, plate_identifier): 
    plate_1 = []           
    for i in data_segregation.columns:
        if plate_identifier in i:
            plate_1.append(i)     
    data_segregation_plate_1 = data_segregation[plate_1]
    return data_segregation_plate_1

# create list of potential wells present
def create_plate_layout():
    alphabet = ['A','B','C','D','E','F','G','H']
    numbers = ['01','02','03','04','05','06','07','08','09','10','11','12']
    plate_annotation = []
    for i in numbers:
        for j in alphabet:
            plate_annotation.append(j+i)
    return plate_annotation

# import segregation table and relabel cols
def import_segregation_table(input_table):
    seg_table = input_table.set_index(['chrom','start','stop'])
    plate_annotation = create_plate_layout()
    seg_table.columns = (list_sample_wells_present(seg_table,plate_annotation))
    seg_table_cols = seg_table.columns
    print(seg_table_cols)
    return seg_table,seg_table_cols

#rename dataframe with unique well indexes more for ease of use
def list_sample_wells_present(data_segregation_plate,plate_annotation):
    
    sample_headers = list(data_segregation_plate.columns)
    for i in plate_annotation:
        for j in  range(len(sample_headers)):
            if i in sample_headers[j]:
                sample_headers[j] = i            
    data_segregation_plate.columns = list(sample_headers)
    return sample_headers

#take a segregation table df and convert to jaccard distance table between samples 
#returns a labeled and indexed df of dist matrix
#neccessary for seaborn
def convert_to_jaccard_dists(seg_table,seg_table_cols):    
    values = []
    for i in range(len(list(seg_table.columns))):
        for j in range(len(list(seg_table.columns))):
            values.append(jaccard(np.array(seg_table.iloc[:,i]),np.array(seg_table.iloc[:,j])))
        
    values = np.array(values).reshape(len(seg_table.columns),len(seg_table.columns))
    values_df = (pd.DataFrame(values))
    values_df.columns = seg_table_cols
    values_df['index'] = seg_table_cols
    values_df.set_index('index',inplace=True)
    return values_df

# process the dist matrix and cluster
#returns sns plotas as a png, the raw and clustered jaccard dist df as a txt
# returns the wells which have over a certain threshold similarity to each other as a txt
def process_and_export_results(values_df,seg_table_cols,similarity_threshold,dpi_value,fig_output_destination,output_raw_df,bad_well_file_destination):
    
    fig = sns.clustermap(values_df)
    fig.savefig(fig_output_destination,dpi=dpi_value)

    new_order = []
    for i in fig.dendrogram_row.reordered_ind:
        new_order.append(seg_table_cols[i])

    clustered_values_df = values_df[new_order]
    clustered_values_df = clustered_values_df.reindex(new_order)
    clustered_values_df.to_csv(output_raw_df,sep = '\t')
    boolean_df = (clustered_values_df > similarity_threshold).astype(int)
    

    for i in range(len(boolean_df)):
        boolean_df.iloc[i,i] = 1

    bad_wells = []
    for i in seg_table_cols:
        if 0 in list(boolean_df.loc[i]):
            bad_wells.append(i)
    print(bad_wells)
    
    with open(bad_well_file_destination, "w") as output:
        output.write(str(bad_wells))

# function wrapper for all functions above
# only input table is neccessary as a str of the filename
# can supply plate identifier if segregation file contains more than one plate
# other parameters are changeable
# currently outputs to CWD

similarity_threshold = 0.2

def contamination_bias_wrapper(input_table, similarity_threshold, folder_to_save_outcome, \
                               plate_identifier = None, dpi_value= 500):
     
    seg_table,seg_table_cols = import_segregation_table(input_table)
    values_df = convert_to_jaccard_dists(seg_table, seg_table_cols)
    process_and_export_results(values_df, seg_table_cols, similarity_threshold, dpi_value,
                      folder_to_save_outcome + plate_identifier + '.figure.png',
                      folder_to_save_outcome + plate_identifier + '.raw_values.txt',
                      folder_to_save_outcome + plate_identifier + '.bad_well.txt')

outcome_folder = 'Plate_cross_contamination/'
import os
if not os.path.exists(outcome_folder):
    os.makedirs(outcome_folder)

contamination_bias_wrapper(GAM_plate_segtable, similarity_threshold, \
                           outcome_folder, 'GAM_plate_segtable_identifier' + '.threshold_' + str(similarity_threshold))

