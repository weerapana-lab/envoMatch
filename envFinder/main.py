
import sys
import os
import argparse
from multiprocessing import Pool
from multiprocessing import cpu_count
import functools

from tqdm import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics.mass import Composition

import modules as src

def parseArgs():
    parser = argparse.ArgumentParser(prog='envFinder',
                                     parents=[src.PARENT_PARSER])

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


def read_pep_stats(fname):
    #cols = ['sequence', 'parent_mz', 'charge', 'scan', 'precursor_scan', 'parent_file', 'sample_name']
    pepStats = pd.read_csv(fname, sep='\t')
    pepStats = pepStats[pepStats['is_modified'].apply(bool)]
    #pepStats = pepStats[cols].drop_duplicates()
    pepStats.reset_index(inplace = True, drop = True)
    # pepStats['precursor_file'] = pepStats['parent_file'].apply(lambda x: x.replace('.ms2', '.ms1'))
    return(pepStats)


def _mkdir(path):
    if not os.path.isdir(path):
        sys.stdout.write('Creating {} ...'.format(path))
        os.mkdir(path)
        sys.stdout.write('Done!\n')


def write_pep_stats(dat, ifname, overwrite = False):
    s = os.path.splitext(ifname)
    base = '{}{}{}'.format(s[0], '_env' if not overwrite else '', s[1])
    sys.stdout.write('\nWriting {}\n'.format(base))
    dat.to_csv(base, sep = '\t', index = False)


def make_of_seq(seq, mod_str = '*'):
    if seq[0] in mod_str:
        raise RuntimeError('{} is an invalid sequence!'.format(seq))

    ret = str()
    for s in seq:
        if s in mod_str:
            ret = ret[:len(ret) - 1] + ret[-1].lower()
        else: ret += s
    return ret


def _annotate_ms1(row, ms1_files=None, args=None, atom_table=None):

    # check that required args are suplied
    assert sum([0 if x is None else 1 for x in [ms1_files, args, atom_table]]) == 3
    ret = False
    _verbose = args.verbose

    #get cit and arg envelopes
    mono_mzs = dict()
    envs = dict()
    charge = row['charge']
    sequence = row['sequence'].upper()
    sequences = [sequence.replace('*', '', x) for x in range(0, sequence.count('*') + 1)]
    for sequences_i, s in enumerate(sequences):
        if args.formula_source == 'calculate':
            comp_temp = atom_table.getComposition(s, row['charge'])
        else:
            mod_diff = atom_table.getComposition(''.join(['*' for _ in range(sequences_i)]), 0, nTerm=False, cTerm=False)
            comp_temp = Composition(formula=src.utils.formula_to_pyteomics_formula(row['formula']), charge=row['charge'])
            comp_temp -= mod_diff

        envs[s] = src.getEnvelope(comp_temp, threshold = 0.01)
        mono_mass = atom_table.getMass(composition=comp_temp, charge=0)
        charge = row['charge']
        mono_mzs[s] = (mono_mass + charge) / charge

    spec = ms1_files[row['parent_file']].getSpectra(row['precursor_scan'],
                                                    (min(mono_mzs.values()) - 5, max(mono_mzs.values()) + 5))

    if spec is None:
        if _verbose:
            sys.stderr.write('Scan: {} not found in {}\n'.format(row['precursor_scan'], row['parent_file']))
        return 'ERROR: Spectrum not found!', 0

    consensus = dict()
    best_score = 0
    best_index = None
    min_mz = sys.float_info.max
    max_mz = 0
    for j, s in enumerate(sequences):
        env = [src.DataPoint(mz, k) for mz, k in envs[s]]
        spec_temp = [src.DataPoint(mz, k) for mz, k in zip(spec['mz'], spec['int'])]

        consensus[s] = src.ConsensusEnvelope(spec_temp, env, sequence = s)
        consensus[s].set_mono(mono_mz = mono_mzs[s])
        consensus[s].annotate(remove_unlabeled=False, verbose=_verbose)

        if consensus[s].envScore > best_score:
            best_score = consensus[s].envScore
            best_index = j
        min_mz = min(min_mz, envs[s][0][0])
        max_mz = max(max_mz, envs[s][-1][0])

    goodEnvTemp = False
    if best_index is not None:
        goodEnvTemp = sequences[best_index] == sequence and best_score >= args.env_co
    ret = goodEnvTemp

    if args.plotEnv:

        fig, ax = plt.subplots(len(sequences), 1, sharex=True)
        fig.set_size_inches(6.4, 2.4 * len(sequences))
        for k, s in enumerate(sequences):
            consensus[s].plotEnv(ax[k], isBest=(k == best_index))

        plt.xlabel('m/z')
        mult = args.mz_step_margin
        plt.xlim(min_mz - (mult / charge) * mult, max_mz + (mult / charge) * mult)

        good_bad_temp = ''
        if args.splitPlots:
            good_bad_temp = 'good' if ret else 'bad'
        path_temp = '{}/envelopes/{}'.format(os.getcwd(), good_bad_temp)
        ofname = '{}/{}_{}_{}_{}.pdf'.format(path_temp,
                                             os.path.splitext(row['parent_file'])[0],
                                             make_of_seq(row['sequence']),
                                             row['scan'],
                                             row['charge'])

        if args.verbose:
            sys.stdout.write('Creating {}...'.format(ofname))

        plt.savefig(ofname)
        plt.close('all')

        if args.verbose:
            sys.stdout.write('Done\n')

    return ret, consensus[sequence].envScore


def main():
    args = parseArgs()

    # calc nThread
    _nThread = args.nThread
    if args.parallel and args.nThread is None:
        _nThread = cpu_count()
    elif not args.parallel and args.nThread is None:
        _nThread=1

    #done once
    pepStats = read_pep_stats(args.ionFinder_output)
    pepStats['good_envelope'] = False

    atom_table = src.AtomTable(args.atom_table)
    if not atom_table.read():
        exit()

    ms1_files = dict()
    ms1_prefix = [os.path.dirname(os.path.abspath(args.ionFinder_output))] + args.ms1_prefix
    for f in pepStats['parent_file'].unique():
        canidate_paths = ['{}/{}.{}'.format(x, os.path.splitext(f)[0], args.file_type) for x in ms1_prefix]
        path = None
        for c in canidate_paths:
            if os.path.isfile(c):
                sys.stdout.write('Found ms1 file for {}!\n\tReading {}'.format(f, c))
                path_temp = c
                break

        if path_temp is None:
            raise RuntimeError('Could not find parent ms1 file for {}!'.format(f))
        else:
            ms1_files[f] = src.Ms1File(path_temp, file_type=args.file_type)
            sys.stdout.write('\n\tDone\n')

    if args.plotEnv:
        path_temp = '{}/envelopes'.format(os.getcwd())
        _mkdir(path_temp)
        if args.splitPlots:
            for s in ['good', 'bad']:
                _mkdir('{}/{}'.format(path_temp, s))

    sys.stdout.write('Searching for envelopes using {} thread(s)...\n'.format(_nThread))
    nRow = len(pepStats.index)
    env_data = list()
    input_lst = [x[1] for x in pepStats.iterrows()]
    show_bar = not(args.verbose and args.parallel == 0)
    if show_bar:
        with Pool(processes=_nThread) as pool:
            env_data = list(tqdm(pool.imap(functools.partial(_annotate_ms1,
                                                                   ms1_files=ms1_files,
                                                                   args=args,
                                                                   atom_table=atom_table),
                                                 input_lst),
                                                 total=len(pepStats.index),
                                                 miniters=1,
                                                 file=sys.stdout))
    else:
        for i, row in pepStats.iterrows():
            sys.stdout.write('\tWorking on {} of {}\r'.format(i, nRow))
            env_data.append(_annotate_ms1(row, ms1_files, args, atom_table))

    pepStats['good_envelope'] = [x[0] for x in env_data]
    pepStats['env_score'] = ['NA' if np.isnan(x[1]) else x[1] for x in env_data]

    sys.stdout.write('\nDone!\n')
    write_pep_stats(pepStats, args.ionFinder_output, overwrite = args.overwrite)


if __name__ == "__main__":
    main()
