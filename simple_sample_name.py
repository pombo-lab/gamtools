import argparse
import os
import sys

parser = argparse.ArgumentParser(description='Takes in any file with paths in the first column and convert them to sample names')
parser.add_argument('input_file', metavar='INPUT_FILE', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='Path to input_file')

def extract_sample_name(path):
    return os.path.basename(path).split('.')[0]

def clean_line(line):
    line = line.strip()
    fields = line.split()
    fields[0] = extract_sample_name(fields[0])
    return '\t'.join(fields)

if __name__ == '__main__':
    args = parser.parse_args()

    for line in args.input_file:
        print clean_line(line)

