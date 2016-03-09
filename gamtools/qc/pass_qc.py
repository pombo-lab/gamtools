import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Print the names of samples passing qc.')
parser.add_argument('stats_file', help='Input statistics file')

args = parser.parse_args()

stats_data = pd.read_csv(args.stats_file, index_col=0, delim_whitespace=True)

enough_mapped_reads = (stats_data.Reads_mapped / stats_data.Reads_sequenced) > 0.15

if 'Type' in stats_data.columns:
    right_type = (stats_data.Type == '1NP')
else:
    sys.stderr.write('No "Type" column found, assuming all samples are NPs\n')
    right_type = True

missing_contam_cols = [col for col in ['Mouse_multiple', 'Mouse_single', 'Human'] if col not in stats_data.columns]
if len(missing_contam_cols) == 0:
    # not_contaminated = ((stats_data.Mouse_multiple + stats_data.Mouse_single) > stats_data.Human)
    not_contaminated = True
else:
    sys.stderr.write('Missing QC columns: {}, check your fastqc_screen configuration\n'.format(missing_contam_cols))
    not_contaminated = True

pd.Series(stats_data[right_type &
           not_contaminated &
           enough_mapped_reads].index).to_csv(sys.stdout, index=False)
