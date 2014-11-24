import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Print the names of samples passing qc.')
parser.add_argument('stats_file', help='Input statistics file')

args = parser.parse_args()

stats_data = pd.read_csv(args.stats_file, index_col=0, delim_whitespace=True)

right_type = (stats_data.Type == '1NP')

enough_mapped_reads = (stats_data.Unique_reads_mapped / stats_data.Reads_sequenced) > 0.05

not_contaminated = ((stats_data.Mouse_multiple + stats_data.Mouse_single) > stats_data.Human)

pd.Series(stats_data[right_type &
           not_contaminated &
           enough_mapped_reads].index).to_csv(sys.stdout, index=False)
