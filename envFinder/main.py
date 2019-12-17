
import sys
import os

import matplotlib.pyplot as plt
import pandas as pd

import modules as src

def read_pep_stats(fname):
    #cols = ['sequence', 'parent_mz', 'charge', 'scan', 'precursor_scan', 'parent_file', 'sample_name']
    pepStats = pd.read_csv(fname, sep='\t')
    pepStats = pepStats[pepStats['is_modified'].apply(bool)]
    #pepStats = pepStats[cols].drop_duplicates()
    pepStats.reset_index(inplace = True)
    # pepStats['precursor_file'] = pepStats['parent_file'].apply(lambda x: x.replace('.ms2', '.ms1'))
    return(pepStats)


def write_pep_stats(dat, ifname, overwrite = False):
    s = os.path.splitext(ifname)
    base = '{}{}{}'.format(s[0], '_env' if not overwrite else '', s[1])
    sys.stdout.write('Writing {}'.format(base))
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


def main():
    args = src.parseArgs()

    #done once
    pepStats = read_pep_stats(args.ionFinder_output)
    pepStats['good_envelope'] = False
    #pepStats = pepStats[pepStats['precursor_scan'] == 8468]
    #pepStats = pepStats[pepStats['precursor_scan'] == 3988]

    atom_table = src.AtomTable(args.atom_table)
    if not atom_table.read():
        exit()

    ms1_files = dict()
    ms1_prefix = [os.path.dirname(os.path.abspath(args.ionFinder_output))] + args.ms1_prefix
    for f in pepStats['parent_file'].unique():
        canidate_paths = ['{}/{}.{}'.format(x, os.path.splitext(f)[0], args.file_type) for x in ms1_prefix]
        path_temp = None
        for c in canidate_paths:
            if os.path.isfile(c):
                sys.stdout.write('Found {}!\n\tReading {}'.format(f, c))
                path_temp = c
                break

        if path_temp is None:
            raise RuntimeError('Could not find parent ms1 file for {}!'.format(f))
        else:
            ms1_files[f] = src.Ms1File(path_temp, file_type=args.file_type)
            sys.stdout.write('\n\tDone\n')

    if args.plotEnv:
        path_temp = '{}/envelopes'.format(os.getcwd())
        if not os.path.isdir(path_temp):
            os.mkdir(path_temp)

    sys.stdout.write('Searching for envelopes...\n')
    nRow = len(pepStats.index)
    for i, row in pepStats.iterrows():
        sys.stdout.write('\tWorking on {} of {}\r'.format(i, nRow))

        #get cit and arg envelopes
        mono_mzs = dict()
        envs = dict()
        charge = row['charge']
        sequence = row['sequence']
        sequences = [sequence.replace('*', '', i) for i in range(0, sequence.count('*') + 1)]
        for s in sequences:
            comp_temp = atom_table.getComposition(s, row['charge'])
            envs[s] = src.getEnvelope(comp_temp, threshold = 0.01)
            mono_mass = atom_table.getMass(s)
            charge = row['charge']
            mono_mzs[s] = (mono_mass + charge) / charge

        spec = ms1_files[row['parent_file']].getSpectra(row['precursor_scan'],
                                                       (min(mono_mzs.values()) - 5, max(mono_mzs.values()) + 5))

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
            consensus[s].annotate(remove_unlabeled=False)

            if consensus[s].envScore > best_score:
                best_score = consensus[s].envScore
                best_index = j
            min_mz = min(min_mz, envs[s][0][0])
            max_mz = max(max_mz, envs[s][-1][0])

        goodEnvTemp = False
        if best_index is not None:
            goodEnvTemp = sequences[best_index] == sequence and best_score >= args.env_co
        #row['good_envelope'] = goodEnvTemp
        pepStats.at[i, 'good_envelope'] = goodEnvTemp

        if args.plotEnv:

            fig, ax = plt.subplots(len(sequences), 1, sharex=True)
            fig.set_size_inches(6.4, 2.4 * len(sequences))
            for i, s in enumerate(sequences):
                consensus[s].plotEnv(ax[i], isBest=(i == best_index))

            plt.xlabel('m/z')
            mult = args.mz_step_margin
            plt.xlim(min_mz - (mult / charge) * mult, max_mz + (mult / charge) * mult)

            ofname = '{}/envelopes/{}_{}_{}_{}.pdf'.format(os.getcwd(),
                                                           row['sample_name'],
                                                           make_of_seq(row['sequence']),
                                                           row['scan'],
                                                           row['charge'])
            plt.savefig(ofname)
            plt.close('all')

    sys.stdout.write('\nDone!\n')
    write_pep_stats(pepStats, args.ionFinder_output, overwrite = args.overwrite)


if __name__ == "__main__":
    main()
