import argparse
import sys
import time
import numpy as np
import itertools
from GamTools import segmentation

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def get_npz_path(args):

    return args.segmentation_file + '.tads_dprime.' + args.region

def tad_name(row):
    return '{0}:{1}-{2}'.format(row.chrom, row.start, row.stop)

def chrom_main(args):

    data = segmentation.open_segmentation(args.segmentation_file)
    tads = segmentation.open_segmentation(args.tad_segmentation_file)

    tad_list = segmentation.region_from_location_string(
        tads, args.region).reset_index()[['chrom', 'start', 'stop']]

    tad_strings = list(tad_list.apply(tad_name, axis=1))

    print 'starting calculation for', args.region
    start_time = time.clock()

    tad_matrix = np.zeros((len(tad_strings),)*3) * np.NaN

    for tad_indexes in itertools.combinations(range(len(tad_strings)),3):

        tad1, tad2, tad3 = tad_indexes
        
        tad_coseg = segmentation.get_matrix_from_regions(
            segmentation.region_from_location_string(data, tad_strings[tad1]),
            segmentation.region_from_location_string(data, tad_strings[tad2]),
            segmentation.region_from_location_string(data, tad_strings[tad3]))
        
        tads_value = np.nanmean(tad_coseg)
        
        for index_comb in itertools.permutations(tad_indexes):
            
            tad_matrix[index_comb] = tads_value

    npz_path = get_npz_path(args)
    np.savez_compressed(npz_path, dprime=tad_matrix, windows=np.array(tad_strings))

    print 'region size is: {0} x {1}'.format(*tad_matrix.shape), 
    print 'Calculation took {0}s'.format(time.clock() - start_time)
    print 'Done!'

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('-r','--region', metavar='REGION', required=True, help='Specific genomic region to calculate matrices for')
parser.add_argument('-s', '--segmentation-file', required=False, help='A segmentation file to use as input')
parser.add_argument('-t', '--tad-segmentation-file', required=False, help='A tad segmentation file to use as input')
parser.set_defaults(func=chrom_main)

if __name__ == "__main__":

    args = parser.parse_args()

    chrom_main(args)

