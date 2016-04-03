import os
import pandas as pd

def screen_out_path(input_fastq):
    
    path_parts = input_fastq.split('.')

    if path_parts[-1] == 'gz':
        path_parts.pop()

    if path_parts[-1] in ['txt','seq','fastq','fq']:
        path_parts.pop()

    path_parts[-1] += '_screen'
    path_parts.append('txt')

    return '.'.join(path_parts)

def get_sample_from_screen_path(fastq_screen_path):
    
    sample = os.path.basename(fastq_screen_path).split('.')[0]
    if sample[-7:] == '_screen':
        sample = sample[:-7]

    return sample

def is_fq_screen_header_row(fields):

    if (len(fields) == 0) or (fields[0][0] == '#') or (fields[0] == 'Library'):
        return True
    else:
        return False

def process_fastqc_screen_line(line):

    fields = line.strip().split()

    if is_fq_screen_header_row(fields):
        row_results = {}

    elif fields[0] == '%Hit_no_libraries:':
        row_results = {'Unmapped': float(fields[1])}
    else:
        row_results = {
            fields[0] + '_single': int(fields[4]) + int(fields[8]),
            fields[0] + '_multiple': int(fields[6]) + int(fields[10]),
            'num_reads': int(fields[1]),
        }

    return row_results

def fq_screen_mapped_other(results_dict):

    mapped_to_other = 0

    for k, reads in results_dict.items():

        parts = k.split('_')
        if (parts[-1] not in ['single', 'multiple']) or (parts[0] in ['Human', 'Mouse']):
            continue

        mapped_to_other += reads

    return mapped_to_other

def get_contamination_stats(fastqc_screen_output_files):

    sample_contamination = []

    for screen_file in fastqc_screen_output_files:

        results = {}

        with open(screen_file) as screen_results:
            for line in screen_results:
                try:
                    results.update(process_fastqc_screen_line(line))
                except ValueError:
                    raise ValueError('Malformed line: "{}"'.format(line))


        total_reads = results['num_reads']

        results['Unmapped'] = int(results['Unmapped'] / 100.0 * total_reads)
        results['Other'] = fq_screen_mapped_other(results)

        for k in results.keys():
            results[k] = 100 * float(results[k]) / total_reads

        results['Sample'] = get_sample_from_screen_path(screen_file)

        sample_contamination.append(results)

    contam_df = pd.DataFrame(sample_contamination)

    contam_df['Human'] = contam_df['Human_single'] + contam_df['Human_multiple']

    return contam_df[['Sample', 'Mouse_single', 'Mouse_multiple', 'Human', 'Other', 'Unmapped']]

def write_contamination_stats(input_files, output_file):

    contam_df = get_contamination_stats(input_files)

    contam_df.to_csv(output_file, sep='\t', index=False)

def contamination_from_doit(dependencies, targets):

    assert len(targets) == 1
    write_contamination_stats(dependencies, targets[0])
