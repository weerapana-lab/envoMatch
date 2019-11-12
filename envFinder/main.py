
import sys
import pdb

import matplotlib.pyplot as plt
import pandas as pd

import modules as src

def read_pep_stats(fnane):
    cols = ['sequence', 'parent_mz', 'charge', 'scan', 'precursor_scan', 'parent_file']
    pepStats = pd.read_csv(fnane, sep='\t')
    pepStats = pepStats[pepStats['is_modified'].apply(bool)]
    pepStats = pepStats[cols].drop_duplicates()
    pepStats.reset_index(inplace = True)
    pepStats['precursor_file'] = pepStats['parent_file'].apply(lambda x: x.replace('.ms2', '.ms1'))
    return(pepStats)


def plotEnv(ax, env):
    ax.stem([x[0] for x in env], [x[1] for x in env], markerfmt=' ')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for x, y in env:
        if y > 0.01:
            ax.annotate('{}'.format(round(x, 4)), xy=(x, y))

    return ax


def main():
    args = src.parseArgs()

    #done once
    pepStats = read_pep_stats('/Volumes/Data/msData/ms2_anotator/citFinder/rorpad_mouse/peptide_cit_stats.tsv')
    #pepStats = pepStats[pepStats['precursor_scan'] == 8468]
    #pepStats = pepStats[pepStats['precursor_scan'] == 5540]

    #done once
    table_keys = ['Citurlline', 'Arginine']
    cit_tableFname = '/Users/Aaron/local/envFinder/db/atom_tables/cit_diff_mod_atoms.txt'
    arg_tableFname = '/Users/Aaron/local/envFinder/db/atom_tables/default_residue_atoms.txt'
    atomTables = dict()
    atomTables[table_keys[0]] = src.AtomTable(cit_tableFname)
    atomTables[table_keys[1]] = src.AtomTable(arg_tableFname)
    for v in atomTables.values():
        if not v.read():
            exit()

    #replace
    # baseDir = '/Volumes/Data/msData/envFinder/ms1/'
    # ms1Dict = dict()
    # for fname in pepStats['precursor_file'].drop_duplicates().to_list():
    #     ms1Dict[fname] = ms1.read('{}/{}'.format(baseDir, fname), use_index=True)

    #ms1File = src.Ms1File('/Volumes/Data/msData/envFinder/ms1/ror_pad_short.ms1')
    ms1File = src.Ms1File('/Volumes/Data/msData/envFinder/ms1/ror_pad.ms1')

    nRow = len(pepStats.index)
    for i, row in pepStats.iterrows():
        #print('Working on {} of {}'.format(i, nRow))

        #get cit and arg envelopes
        mono_mzs = dict()
        envs = dict()
        for k, v in atomTables.items():
            comp_temp = v.getComposition(row['sequence'], row['charge'])
            envs[k] = src.getEnvelope(comp_temp, threshold = 0.01)
            mono_mass = v.getMass(row['sequence'])
            charge = row['charge']
            mono_mzs[k] = (mono_mass + charge) / charge

        spec = ms1File.getSpectra(row['precursor_scan'],
                                  (min(mono_mzs.values()) - 5, max(mono_mzs.values()) + 5))

        #spec = ms1File.getSpectra(row['precursor_scan'], Ms1File.getMzRange())

        #actualEnv = src.ms1.findEnvelope(env, spec)

        # cEnv = consensusEnvelope.ConsensusEnvelope()
        # cEnv.initalize(env, 'cit')
        #
        # cEnv.add(actualEnv)

        charge = row['charge']
        consensus = dict()
        best_score = 0
        best_index = None
        min_mz = sys.float_info.max
        max_mz = 0
        for i, k in enumerate(table_keys):
            env = [src.DataPoint(mz, i) for mz, i in envs[k]]
            spec_temp = [src.DataPoint(mz, i) for mz, i in zip(spec['mz'], spec['int'])]

            consensus[k] = src.ConsensusEnvelope(spec_temp, env, sequence = row['sequence'])
            consensus[k].set_mono(mono_mz = mono_mzs[k])
            consensus[k].annotate(remove_unlabeled=False)

            if consensus[k].envScore > best_score:
                best_score = consensus[k].envScore
                best_index = i
            min_mz = min(min_mz, envs[k][0][0])
            max_mz = max(max_mz, envs[k][-1][0])

        fig, ax = plt.subplots(len(table_keys), 1, sharex = True)
        for i, k in enumerate(table_keys):
            consensus[k].plotEnv(ax[i], isBest = (i == best_index))

        plt.xlabel('m/z')
        mult = args.mz_step_margin
        plt.xlim(min_mz - (mult / charge) * mult, max_mz + (mult / charge) * mult)
        #plt.show()
        print('{}\t{}'.format(row['sequence'], atomTables[table_keys[0]].getMass(row['sequence'])))

        #fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
        #ax1 = plotEnv(ax1, env)
        #ax2 = plotEnv(ax2, actualEnv)

        plt.savefig('/Users/Aaron/local/envFinder/testFiles/env_test/test_env_{}.pdf'.format(row['precursor_scan']))
        plt.close('all')



if __name__ == "__main__":
    main()
