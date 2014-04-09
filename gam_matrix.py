from GamTools import GamExperiment
import argparse
import os
import sys
import time

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('segmentation_file', help='A segmentation file to use as input')
parser.add_argument('-p','--num-processes', metavar='NUM_PROC', type=int, default=1, help='Number of processes to use to calculate matrix')
parser.add_argument('-z','--compression', metavar='COMPRESSION', help='Level of compression to use in hdf5 files', default=None)
parser.add_argument('-c','--chromosomes', metavar='CHROMOSOME', nargs='+', help='Specific chromosomes to calculate matrices for (default is all chromosomes not ending in _random)')
parser.add_argument('-o','--overwrite', action='store_true', help='Specific chromosomes to calculate matrices for (default is all chromosomes not ending in _random)')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-s','--hdf5-path', metavar='hdf5_path', help='path to store matrix in hdf5 format')
group.add_argument('-a','--auto-hdf5', action='store_true',
                    help='automatically determine path to store matrix in hdf5 format')

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

    base = os.path.basename(segmentation_path)
    return base + '.matrix.hdf5'

def main(args):

    if args.auto_hdf5:
        args.hdf5_path = get_auto_hdf5_path(args.segmentation_file)

    if os.path.exists(args.hdf5_path):
        if args.overwrite or query_yes_no('This hdf5 file already exists, are you sure you want to erase it?'):
            os.unlink(args.hdf5_path)
        else:
            sys.exit('Please specify a path to create a new hdf5 file')

    exp = GamExperiment.from_multibam(args.segmentation_file, args.hdf5_path, args.num_processes,
                                     compression=args.compression)

    if not args.chromosomes:
        args.chromosomes = [ c for c in sorted(list(set(map(lambda i: i[0], list(exp.experimental_data.windows))))) if not c[-7:] == '_random' ]

    print 'using {0} processes'.format(args.num_processes)
    for chrom in args.chromosomes:
        print 'starting calculation for', chrom
        start_time = time.clock()
        chrom_shape = exp.frequencies(chrom).shape
        print 'chrom size is: {0} x {1}'.format(*chrom_shape), 
        print 'Calculation took {0}s'.format(time.clock() - start_time)
    print 'Done!'

if __name__ == "__main__":

    args = parser.parse_args()

    main(args)
