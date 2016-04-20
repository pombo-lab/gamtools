from . import segregation

def select_samples(segregation_path, sample_names, output_file, drop=False):

    segregation_data = segregation.open_segregation(segregation_path)

    column_mapping = segregation.map_sample_name_to_column(segregation_data)
    column_names = [column_mapping[name] for name in sample_names]

    if drop:
        subset = segregation_data.drop(column_names, axis=1)
    else:
        subset = segregation_data[column_names]

    subset.to_csv(output_file, index=True, sep='\t')

def select_samples_from_file(segregation_path, sample_names_path, output_file, drop=False):

    names = []

    with open(sample_names_path, 'r') as sample_names_file:
        for line in sample_names_file:
            names.append(line.strip())

    select_samples(segregation_path, names, output_file, drop)

def select_samples_from_args(args):

    select_samples(args.segregation_file, args.sample_names, args.output_file, args.drop_samples)

def select_samples_from_doit(dependencies, targets):

    for dep_file in dependencies:
        dep_ext = dep_file.split('.')[-1]

        if dep_ext == 'multibam':
            segregation_path = dep_file
        elif dep_ext == 'txt':
            sample_names_path = dep_file

    assert len(targets) == 1

    select_samples_from_file(segregation_path, sample_names_path, list(targets)[0])
