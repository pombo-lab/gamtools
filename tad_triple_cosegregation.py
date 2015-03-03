import argparse
import sys
import time
import numpy as np
import itertools
from GamTools import segmentation
from GamTools.coseg_3 import coseg_all_3

def get_npz_path(args):

    return args.segmentation_file + '.tads_dprime.' + '^'.join(args.regions)

def tad_name(row):
    return '{0}:{1}-{2}'.format(row.chrom, row.start, row.stop)

def chrom_main(args):

    data = segmentation.open_segmentation(args.segmentation_file)
    tads = segmentation.open_segmentation(args.tad_segmentation_file)

    def get_tads(loc):
        
        chrom = loc.split(':')[0]
        
        chrom_df = segmentation.region_from_location_string(tads, chrom).reset_index()[['chrom', 'start', 'stop']].reset_index()
        chrom_df['name'] = chrom_df.apply(tad_name, axis=1)
        chrom_df = chrom_df.set_index(['chrom', 'start', 'stop'])
        
        tad_strings = segmentation.region_from_location_string(chrom_df, loc)
        
        return tad_strings.itertuples(index=False)

    def get_region_tad_combs(loc1, loc2, loc3):
            
        tads1, tads2, tads3 = [get_tads(loc) for loc in 
                    (loc1, loc2, loc3)]

        tad_combs = itertools.product(tads1, tads2, tads3)

        uniq_combs = list(set([tuple(sorted(comb)) for comb in tad_combs]))

        return uniq_combs

    tad_combinations = get_region_tad_combs(*args.regions)

    print 'starting calculation for', '^'.join(args.regions)
    start_time = time.clock()

    result = []
    res_append = result.append

    for tad_indexes, tad_strings in [zip(*comb) for comb in tad_combinations]:

        if len(set(tad_indexes)) != 3:
            continue

        tad_data = [np.array(segmentation.region_from_location_string(data, tad_str)).astype(int) for
                    tad_str in tad_strings]

        tad_coseg = coseg_all_3(*tad_data)

        tads_value = np.nanmean(tad_coseg)
        
        res_append(list(tad_indexes) + [tads_value])

    npz_path = get_npz_path(args)
    np.savez_compressed(npz_path, dprime=result)

    print 'region size is: {0}'.format(len(result)), 
    print 'Calculation took {0}s'.format(time.clock() - start_time)
    print 'Done!'

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-r','--regions', metavar='REGION', nargs='+', required=True, help='Specific genomic regions to calculate matrices for')
parser.add_argument('-s', '--segmentation-file', required=False, help='A segmentation file to use as input')
parser.add_argument('-t', '--tad-segmentation-file', required=False, help='A tad segmentation file to use as input')
parser.set_defaults(func=chrom_main)

if __name__ == "__main__":

    args = parser.parse_args()

    assert len(args.regions) == 3

    chrom_main(args)

