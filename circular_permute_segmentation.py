import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Circularly permute every column of a segmentation file')
parser.add_argument('segmentation_file', help='A segmentation file to use as input')

def read_data(fpointer):
    windows = []
    data = []
    for lnum, line in enumerate(fpointer):
        if lnum == 0:
            header = line.strip()
        else:
            fields = line.strip().split()
            windows.append(fields[:2])
            data.append(map(int,fields[2:]))

    data = np.array(data).transpose()

    return header, data, windows

def permute_column(column):
    no_windows = len(data[0])
    offset = np.random.randint(no_windows)
    permutation = list(data[0][offset:])
    permutation.extend(data[0][:offset])
    return permutation

def permute_data(data):
    new_array = []
    for i in range(len(data)):
        new_array.append(permute_column(data[i]))
    return np.array(new_array)

def output_permutation(header, permutation, windows):
    print header
    for i,row in enumerate(permutation.transpose()):
        full_row = windows[i] + map(str, list(row))
        print '\t'.join(full_row)

if __name__ == '__main__':

    args = parser.parse_args()
    
    with open(args.segmentation_file) as fpointer:
        header, data, windows = read_data(fpointer)

    permutation = permute_data(data)

    output_permutation(header, permutation, windows)
