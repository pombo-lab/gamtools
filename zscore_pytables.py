import tables
import pandas as pd
import itertools
import re
import sys
import numpy as np

class Interaction(tables.IsDescription):
    chrom      = tables.StringCol(2)   # 1-character String
    grid_i    = tables.Int32Col()      # 32-bit integer
    grid_j    = tables.Int32Col()      # 32-bit integer
    score = tables.Float32Col()    

h5file = tables.openFile("zscores.h5", mode = "w", title = "Test file")

group = h5file.createGroup("/", 'interaction', 'Interaction information')

table = h5file.createTable(group, 'zscore', Interaction, "Z scores")

interaction = table.row


def get_mean_std(chrom_data):
    indices = itertools.product(xrange(chrom_data.shape[0]), xrange(chrom_data.shape[1]))
    table_data = [ (p, q, chrom_data[p,q]) for p,q in indices ] 
    df = pd.DataFrame.from_records(table_data, columns=['x', 'y', 'correlation'])
    df['dist'] = np.abs(df.x - df.y)
    dists = df.groupby('dist')['correlation']
    means = dists.mean()
    stdevs = dists.std()
    return means, stdevs

def get_chr_string(npz_path):
    return re.search('chr[0-9X]{1,2}', npz_path).group(0)[3:]

def get_z_scores(npz_path):
    
    chrom_data = np.load(npz_path)['corr']
    
    chrom_name = get_chr_string(npz_path)
    
    means, stdevs = get_mean_std(chrom_data)
    
    def get_z_score(x, dist):
        return (x - means[dist]) / stdevs[dist]
    
    for i in xrange(chrom_data.shape[0]):
        for j in xrange(chrom_data.shape[1]):
            current_z = get_z_score(chrom_data[i,j], np.abs(i-j))
            if np.isfinite(current_z):
                interaction['chrom']  = chrom_name
                interaction['grid_i'] = i
                interaction['grid_j'] = j
                interaction['score'] = current_z
                # Insert a new particle record
                interaction.append()
        table.flush()
    table.flush()


for p in sys.argv[1:]:
    print p
    get_z_scores(p)

h5file.root.interaction.zscore.cols.score.createCSIndex()

table = h5file.root.interaction.zscore

print [ (x['score'], x['grid_i'], x['grid_j']) for x in table.itersorted('score',start=None,stop=15)]

h5file.close()

