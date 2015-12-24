import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('--doit-num-process', metavar='NUM_PROCESSES', 
                    type=int, default=1, 
                    help='Number of tasks to execute in parallel')

parser.add_argument('--doit-db-file', metavar='DB_FILE',
                    default=os.path.join(os.getcwd(), '.doit.db'), 
                    help='Location of doit database file')

parser.add_argument('-v', '--verbosity', action='count',
                    help='Set the verbosity level (-v is least verbose, -vvv is most)')

def config_from_args(args):
    config = {}

    if hasattr(args, 'doit_num_process') and args.doit_num_process is not None:

        config['num_process'] = args.doit_num_process

    if hasattr(args, 'doit_db_file') and args.doit_db_file is not None:

        config['dep_file'] = args.doit_db_file

    if hasattr(args, 'verbosity') and args.verbosity is not None:

        if args.verbosity > 3:
            args.verbosity = 3

        config['verbosity'] = args.verbosity

    return config
