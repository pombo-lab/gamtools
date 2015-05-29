import argparse

parser = argparse.ArgumentParser(description='Normalise a bedgraph file by the number of lines in a bed file')
parser.add_argument('bedgraph', metavar='BEDGRAPH', 
				   help='The bedgraph to be normalised')
parser.add_argument('reads', metavar='NUMBER_OF_READS', 
				   help='The total number of reads in the file')

args = parser.parse_args()

def normalise_line(line):
    fields = line.split()
    fields[3] = str((float(fields[3])/float(args.reads)) * 1000000.0)
    return '\t'.join(fields)

for line in open(args.bedgraph):
    print normalise_line(line)
