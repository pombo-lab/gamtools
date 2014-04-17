from GamTools import GamExperiment
import argparse
import os
import sys
import time
import numpy as np

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

def get_auto_hdf5_path(segmentation_path):

    return segmentation_path + '.matrix.hdf5'

def create_main(args):

    if args.auto_hdf5:
        args.hdf5_path = get_auto_hdf5_path(args.segmentation_file)

    if os.path.exists(args.hdf5_path):
        if args.overwrite or query_yes_no('This hdf5 file already exists, are you sure you want to erase it?'):
            os.unlink(args.hdf5_path)
        else:
            sys.exit('Please specify a path to create a new hdf5 file')

    exp = GamExperiment.from_multibam(args.segmentation_file, args.hdf5_path,
                                     compression=args.compression)

    exp.close()

def chrom_main(args):

    print 'using {0} processes'.format(args.num_processes)
    exp = GamExperiment.from_multibam_no_cache(args.segmentation_file, args.num_processes)
    print 'starting calculation for', args.chromosome
    start_time = time.clock()
    chrom_path = args.segmentation_file + '.chrom.' + args.chromosome
    chrom_freq = exp.frequencies(args.chromosome)
    np.savez_compressed(chrom_path, freqs=chrom_freq)
    chrom_shape = chrom_freq.shape
    print 'chrom size is: {0} x {1}'.format(*chrom_shape), 
    print 'Calculation took {0}s'.format(time.clock() - start_time)
    exp.close()
    print 'Done!'

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
subparsers = parser.add_subparsers()

get_chrom = subparsers.add_parser('chrom')
get_chrom.add_argument('-p','--num-processes', metavar='NUM_PROC', type=int, default=1, help='Number of processes to use to calculate matrix')
get_chrom.add_argument('-c','--chromosome', metavar='CHROMOSOME', required=True, help='Specific chromosome to calculate matrices for')
get_chrom.add_argument('segmentation_file', help='A segmentation file to use as input')
get_chrom.set_defaults(func=chrom_main)

create_options = subparsers.add_parser('create')
create_options.add_argument('-z','--compression', metavar='COMPRESSION', help='Level of compression to use in hdf5 files', default=None)
create_options.add_argument('-o','--overwrite', action='store_true', help='Specific chromosomes to calculate matrices for (default is all chromosomes not ending in _random)')
group = create_options.add_mutually_exclusive_group(required=True)
group.add_argument('-s','--hdf5-path', metavar='hdf5_path', help='path to store matrix in hdf5 format')
group.add_argument('-a','--auto-hdf5', action='store_true',
                    help='automatically determine path to store matrix in hdf5 format')
create_options.set_defaults(func=create_main)

if __name__ == "__main__":

    args = parser.parse_args()

    args.func(args)
