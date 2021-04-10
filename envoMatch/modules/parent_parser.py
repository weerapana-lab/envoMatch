
import os
import argparse

ATOM_TABLE_PATH = os.path.abspath('{}/../../db/atom_tables/cit_diff_mod_atoms.txt'.format(os.path.dirname(os.path.abspath(__file__))))

PARENT_PARSER = argparse.ArgumentParser(add_help=False)

PARENT_PARSER.add_argument('--env_co', default=0.8, type=float,
                           help='Envelope correlation score cutoff.')

PARENT_PARSER.add_argument('--mz_step_margin', default=2, type=int,
                           help='Margin above and below envelope in plot.')

PARENT_PARSER.add_argument('-t', '--file_type', choices=['ms1', 'mzXML', 'mzML'], default='mzXML',
                           help='MS-1 input file type. Default is mzXML.')

PARENT_PARSER.add_argument('--ms1_prefix', action='append',
                           help='Append directory to search path for ms1 files. '
                                'By default only the current working directory is used. '
                                'Additional directories are searched in the order they are provided.')

PARENT_PARSER.add_argument('-f', '--formula_source', choices=['input', 'calculate'], default='calculate',
                           help='Where should peptide formulas come from? Default is calculate.')

PARENT_PARSER.add_argument('-s', '--pre_scan_src', choices=['input', 'ms1'], default='ms1',
                           help='Where should precursor scans come from. '
                           'Chose either the "precursor_scan" column (input) or build precursor list '
                           'from input MS-1 files (ms1). Default is ms1.')

PARENT_PARSER.add_argument('-a', '--atom_table', default=ATOM_TABLE_PATH,
                           help='Path to atom table to use in calculating envelopes. '
                                'Default is: {}'.format(ATOM_TABLE_PATH))

PARENT_PARSER.add_argument('--plotEnv', action='store_true', default=False,
                           help='Should plot of envelopes be saved?')

PARENT_PARSER.add_argument('--splitPlots', action='store_true', default=False,
                           help='Split "good" and "bad" envelope plots in separate directories.')

PARENT_PARSER.add_argument('--parallel', choices=[0, 1], type=int, default=1,
                           help='Chose whether envelope matching should be performed in parallel.'
                                ' Parallel processing is performed on up to the number of logical cores on your system. '
                                '1 is the default.')

PARENT_PARSER.add_argument('--nThread', type=int, default=None,
                           help='Chose how many threads to use for parllel processing. '
                                'This option overrides the --parallel option.')

PARENT_PARSER.add_argument('--overwrite', type=int, choices=[0, 1], default=0,
                           help='Should ionFinder_output be overwritten?')

PARENT_PARSER.add_argument('-v', '--verbose', action='store_true', default=False,
                           help='Print verbose output?')

PARENT_PARSER.add_argument('input_file', type=str,
                           help='ionFinder output file to read.')

