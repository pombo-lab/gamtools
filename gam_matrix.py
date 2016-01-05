import argparse
import sys
import time
import numpy as np
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

    return args.segmentation_file + '.dprime.' + args.region

def chrom_main(args):

    data = segmentation.open_segmentation(args.segmentation_file)
    print 'starting calculation for', args.region
    start_time = time.clock()
    region = segmentation.region_from_location_string(data, args.region)
    matrix = segmentation.get_dprime_from_regions(region)
    windows = np.array(list(region.index))
    npz_path = get_npz_path(args)
    np.savez_compressed(npz_path, scores=matrix, windows=windows)
    region_shape = matrix.shape
    print 'region size is: {0} x {1}'.format(*region_shape), 
    print 'Calculation took {0}s'.format(time.clock() - start_time)
    print 'Done!'

parser = argparse.ArgumentParser(description='Generate GAM heat maps from a table of positive windows by NP.')
parser.add_argument('-r','--region', metavar='REGION', required=True, help='Specific genomic region to calculate matrices for')
parser.add_argument('-s', '--segmentation_file', required=False, help='A segmentation file to use as input')
parser.set_defaults(func=chrom_main)

if __name__ == "__main__":

    args = parser.parse_args()

    chrom_main(args)
