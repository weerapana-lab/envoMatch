
from pyteomics import ms1 as Ms1File
import sys
import matplotlib.pyplot as plt

import pandas as pd

import atom_table
import parse_args
import ms1


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
    args = parse_args.parseArgs()

    #done once
    pepStats = read_pep_stats('/Volumes/Data/msData/ms2_anotator/citFinder/rorpad_mouse/peptide_cit_stats.tsv')
    pepStats = pepStats[pepStats['precursor_scan'] == 8468]

    #done once
    atomTableFname = '/Users/Aaron/local/envFinder/db/atom_tables/cit_diff_mod_atoms.txt'
    atomTable = atom_table.AtomTable(atomTableFname)
    if not atomTable.read():
        exit()

    #replace
    # baseDir = '/Volumes/Data/msData/envFinder/ms1/'
    # ms1Dict = dict()
    # for fname in pepStats['precursor_file'].drop_duplicates().to_list():
    #     ms1Dict[fname] = ms1.read('{}/{}'.format(baseDir, fname), use_index=True)

    ms1File = Ms1File.read('/Volumes/Data/msData/envFinder/ms1/ror_pad_short.ms1', use_index=True)
    #ms1File = Ms1File.read('/Volumes/Data/msData/envFinder/ms1/ror_pad.ms1', use_index=True)

    nRow = len(pepStats.index)
    for i, row in pepStats.iterrows():
        print('Working on {} of {}'.format(i, nRow))

        c = atomTable.getComposition(row['sequence'], row['charge'])
        env = atom_table.getEnvelope(c)


        spec = ms1File.get_by_id(str(row['precursor_scan']).zfill(6))

        

        actualEnv = ms1.findEnvelope(env, spec)

        #fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
        #ax1 = plotEnv(ax1, env)
        #ax2 = plotEnv(ax2, actualEnv)

        #plt.savefig('/Users/Aaron/local/envFinder/testFiles/test_env_{}.pdf'.format(row['scan']))
        #plt.close('all')



if __name__ == "__main__":
    #main()
    drawMassDefectExample()
