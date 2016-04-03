import pandas as pd


def get_reference(condition_ref_string, stats_df):
    if condition_ref_string in stats_df.columns:
        return stats_df.loc[:, condition_ref_string]
    
    try:
        return float(condition_ref_string)
    except ValueError:
        return condition_ref_string


def get_references(left_str, right_str, stats_df):
    left = get_reference(left_str, stats_df)
    right = get_reference(right_str, stats_df)
    
    if not (isinstance(left, pd.Series) or isinstance(right, pd.Series)):
        raise QcParamError(
            """Neither {} nor {} matches a column in the stats file {}.
            Current options are: \n\t'{}'""".format(
                left_str, right_str, stats_df_path,
                "'\n\t'".join(stats_df.columns)))
        
    return left, right


def comparison_from_operator(op, left, right):
    
    if op in ['=', '==', 'eq', 'equals']:
        comparison = (left == right)
    elif op in ['>', 'gt', 'greater_than']:
        comparison = (left > right)
    elif op in ['>=', 'gte', 'greater_than_or_equal_to']:
        comparison = (left >= right)
    elif op in ['<', 'lt', 'less_than']:
        comparison = (left < right)
    elif op in ['<=', 'lte', 'less_than_or_equal_to']:
        comparison = (left <= right)
    elif op in ['!=', 'neq', 'not_equal_to']:
        comparison = (left != right)
    else:
        raise QcParamError(
            'Operator {} not recognized'.format(op))
        
    return comparison


class QcParamError(Exception):
    pass


def do_comparison(left_str, op, right_str, stats_df):
    
    left, right = get_references(left_str, right_str, stats_df)
    
    comparison = comparison_from_operator(op, left, right)
        
    if not isinstance(comparison, pd.Series):
        raise QcParamError('Comparison did not return a series object')
        
    return comparison


def parse_conditions_file(conditions_file, stats_df):

    conditions = []

    for line_no, line in enumerate(conditions_file, 1):
        
        fields = line.split()
        
        if (line[0] == '#') or len(fields) == 0: continue
        
        left_str, op, right_str = fields
        
        try:
            this_comparison =  do_comparison(left_str, op, right_str, stats_df)
        except QcParamError, err_msg:
            raise QcParamError(
                'Error in QC parameters file {}\nLine {}: {}'.format(
                    conditions_file_path, line_no, err_msg))
        
        conditions.append(this_comparison)
    
    return conditions


def samples_passing_qc(conditions_file_path, stats_df_path):

    stats_df = pd.read_csv(stats_df_path, delim_whitespace=True)

    with open(conditions_file_path, 'r') as conditions_file:
        conditions = parse_conditions_file(conditions_file, stats_df)

    sample_passes_qc = pd.concat(conditions, axis=1).all(axis=1)

    return stats_df.loc[sample_passes_qc, 'Sample']

def create_passqc_file(conditions_file_path, stats_df_path, passqc_path):

    samples_passing = samples_passing_qc(conditions_file_path, stats_df_path)
    samples_passing.to_csv(passqc_path, sep='\t', index=False, header=False)

def samples_passing_qc_from_doit(targets, dependencies):

    for file_path in dependencies:

        file_ext = file_path.split('.')[-1]
        if file_ext == 'cfg':
            conditions_file_path = file_path
        elif file_ext == 'txt':
            stats_df_path = file_path

    create_passqc_file(conditions_file_path, stats_df_path, list(targets)[0])
