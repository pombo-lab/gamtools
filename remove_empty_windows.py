import pandas as pd
import sys
data = pd.read_csv(sys.argv[1],
                   delim_whitespace=True, index_col=[0,1], skiprows=1,header=None)
for i,window in enumerate(data.index):
    if any(list(data.iloc[i,:])):
        print window[0] + '-' + window[1]
