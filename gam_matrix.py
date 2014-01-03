from GamTools import GamExperiment
import argparse

parser = argparse.ArgumentParser(description='Calculate coverage over different window sizes for a list of bam files.')
parser.add_argument('segmentation_file', help='A segmentation file to use as input')
parser.add_argument('-s','--hdf5-path', metavar='HDF5_PATH', required=True, help='Path to store matrix in hdf5 format')

def main(args):

    exp = GamExperiment.from_multibam(args.segmentation_file, args.hdf5_path)

    for chrom in sorted(list(set(map(lambda i: i[0], list(exp.experimental_data.data.columns))))):
        if not chrom[-7:] == '_random':
            print chrom
            print exp.get_chrom_processed_matrix(chrom).shape
    print 'Done!'

if __name__ == "__main__":

    args = parser.parse_args()

    main(args)
