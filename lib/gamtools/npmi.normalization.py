def NPMI(A,B):
    import numpy as np
    window_A = sum(A)/len(A)
    window_B = sum(B)/len(B)
    joint_AB = (np.count_nonzero((A+B)> 1))/len(A)
    pmi = np.log2(joint_AB/(window_A * window_B))
    npmi = pmi/-np.log2(joint_AB)
    return(npmi)

def npmi_normalize_matrix(segregation_matrix):
    import numpy as np
    np.warnings.filterwarnings('ignore')
    scores = []
    segregation_matrix = np.array(segregation_matrix)
    for window_A in segregation_matrix:
        for window_B in segregation_matrix:
            scores.append(NPMI(window_A,window_B))
    scores = np.array(scores)
    matrix = scores.reshape(len(segregation_matrix),len(segregation_matrix))
    return pd.DataFrame(matrix)

def get_npmi_scale (matrix):
    import numpy as np
    matrix_con = np.concatenate (matrix)
    matrix_remove_nans = matrix_con[~np.isnan(matrix_con)]
    minscale = 0
    maxscale = np.percentile (matrix_remove_nans, 99)
    return minscale, maxscale

from gamtools import segregation
region = 'chr:start-stop'

npmi_matrix = npmi_normalize_matrix (segregation.region_from_location_string (segregation_table, region))

