import argparse
import os

parser = argparse.ArgumentParser(description='Takes in a sorted bed file and returns the proportion of gaps between segments equal to the smallest gap')
parser.add_argument('segmentation_beds', metavar='SEGMENTATION_BED', nargs='*', help='Path to one or more segmentation bed files (must be sorted)')

def count_gaps(bed_handle):
    gaps = []
    last_end = 0
    last_chrom = False
    for line in bed_handle:
        chrm, start, end = line.strip().split()
        start, end = int(start), int(end)
        if last_chrom == chrm:
            gaps.append(start - last_end)
        last_end = end
        last_chrom = chrm
    return gaps

def proportion_smallest(gaps):
    no_smallest = float(len([gap for gap in gaps if gap == min(gaps)]))
    return no_smallest / len(gaps)

def sample_name(sample_path):
    sample_basename = os.path.basename(sample_path)
    return sample_basename.split('.')[0]

if __name__ == '__main__':
    args = parser.parse_args()

    for bed_path in args.segmentation_beds:
        with open(bed_path, 'r') as bed_handle:
            print sample_name(bed_path), proportion_smallest(count_gaps(bed_handle))
