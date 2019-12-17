
import argparse
import os

ATOM_TABLE_PATH = os.path.abspath('../db/atom_tables/cit_diff_mod_atoms.txt')


def parseArgs():
    parser = argparse.ArgumentParser(prog='envFinder')

    parser.add_argument('--env_co', default = 0.8, type = float,
                        help = 'Envelope correlation score cutoff.')

    parser.add_argument('--mz_step_margin', default = 2, type= int,
                        help = 'Margin above and below envelope in plot.')

    parser.add_argument('-t', '--file_type', choices=['ms1', 'mzXML', 'mzML'], default = 'mzXML',
                        help = 'MS-1 input file type. Default is mzXML.')

    parser.add_argument('-p', '--ms1_prefix', action = 'append',
                        help = 'Append directory to search path for ms1 files. '
                               'By default only the current working directory is used. '
                               'Additional directories are searched in the order they are provided.')

    parser.add_argument('-a', '--atom_table', default = ATOM_TABLE_PATH,
                        help = 'Path to atom table to use in calculating envelopes. '
                               'Default is: {}'.format(ATOM_TABLE_PATH))

    parser.add_argument('--plotEnv', action = 'store_true', default = False,
                        help = 'Should plot of envelopes be saved?')

    parser.add_argument('ionFinder_output', type = str,
                        help='ionFinder output file to read.')


    args = parser.parse_args()

    #check ms1 prefix list
    if args.ms1_prefix is None:
        args.ms1_prefix = [os.path.dirname(os.path.abspath(args.ionFinder_output))]
    else:
        args.ms1_prefix.insert(0, os.path.dirname(os.path.abspath(args.ionFinder_output)))

    for p in args.ms1_prefix:
        if not os.path.isdir(p):
            raise RuntimeError('{} is not a directory!'.format(p))

    return args
