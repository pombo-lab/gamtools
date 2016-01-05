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

    return args.segmentation_file + '.freqs.' + args.name

def chrom_main(args):

    data = segmentation.open_segmentation(args.segmentation_file)
    if args.region is not None:
        region1 = segmentation.region_from_location_string(data, args.region)
        region2 = region1
    else:
        i1_start, i1_stop = args.index1.split('-')
        i2_start, i2_stop = args.index2.split('-')
        i1_start, i1_stop, i2_start, i2_stop = int(i1_start), int(i1_stop), int(i2_start), int(i2_stop)
        region1 = data.iloc[i1_start:i1_stop]
        region2 = data.iloc[i2_start:i2_stop]

    reg1_parts = region1.iloc[0].name + region1.iloc[-1].name
    reg2_parts = region2.iloc[0].name + region2.iloc[-1].name

    print 'starting calculation for {}:{}-{} against {}:{}-{}'.format(reg1_parts[0], reg1_parts[1], reg1_parts[-1],
                                                                      reg2_parts[0], reg2_parts[1], reg2_parts[-1])

    start_time = time.clock()
    freqs = segmentation.get_cosegregation_freqs(region1, region2)
    windows1 = np.array(list(region1.index))
    windows2 = np.array(list(region2.index))
    npz_path = get_npz_path(args)
    np.savez_compressed(npz_path, freqs=freqs, windows1=windows1, windows2=windows2)
    region_shape = freqs.shape
    print 'region size is: {0} x {1}'.format(*region_shape), 
    print 'Calculation took {0}s'.format(time.clock() - start_time)
    print 'Done!'

parser = argparse.ArgumentParser(description='Generate GAM co-segregation heat maps from a table of positive windows by NP.')
parser.add_argument('-r','--region', metavar='REGION', help='Specific genomic region to calculate matrices for', default=None)
parser.add_argument('-i1','--index1', metavar='START-STOP', help='First set of indices to calculate co-segregation for.', default=None)
parser.add_argument('-i2','--index2', metavar='START-STOP', help='Second set of indices to calculate co-segregation for.', default=None)
parser.add_argument('-s', '--segmentation_file', required=True, help='A segmentation file to use as input')
parser.set_defaults(func=chrom_main, name='')

if __name__ == "__main__":

    args = parser.parse_args()
    if args.region is None:
        args.name = '{}_{}'.format(args.index1, args.index2)

    else:
        args.name = args.region

    chrom_main(args)
