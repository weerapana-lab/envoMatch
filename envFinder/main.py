
from pyteomics import ms1
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import ggplot
import pandas as pd

import peptide
import parse_args


def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]


def read_pep_stats(fnane):
    cols = ['sequence', 'parent_mz', 'charge', 'scan', 'precursor_scan', 'parent_file']
    pepStats = pd.read_csv(fnane, sep='\t')
    pepStats = pepStats[pepStats['is_modified'].apply(bool)]
    pepStats = pepStats[cols].drop_duplicates()
    pepStats['precursor_file'] = pepStats['parent_file'].apply(lambda x: x.replace('.ms2', '.ms1'))
    return(pepStats)


def main():
    args = parse_args.parseArgs()

    #done once
    pepStats = read_pep_stats('/Volumes/Data/msData/ms2_anotator/citFinder/rorpad_mouse/peptide_cit_stats.tsv')
    pepStats = pepStats[pepStats['precursor_scan'] == 8468]

    #done once
    atomTableFname = '/Users/Aaron/local/envFinder/db/atom_tables/cit_diff_mod_atoms.txt'
    atomTable = peptide.AtomTable(atomTableFname)
    if not atomTable.read():
        exit()

    #replace
    # baseDir = '/Volumes/Data/msData/envFinder/ms1/'
    # ms1Dict = dict()
    # for fname in pepStats['precursor_file'].drop_duplicates().to_list():
    #     ms1Dict[fname] = ms1.read('{}/{}'.format(baseDir, fname), use_index=True)

    ms1File = ms1.read('/Volumes/Data/msData/envFinder/ms1/ror_pad_short.ms1', use_index=True)

    for i, row in pepStats.iterrows():
        c = atomTable.getComposition(row['sequence'], row['charge'])
        env = peptide.getEnvelope(c)


        spec = ms1File.get_by_id(str(row['precursor_scan']).zfill(6))

        parentMass = 661.33441

        val = find_nearest(spec['m/z array'], parentMass)
        print(abs(parentMass - val))

        fig, ax = plt.subplots()
        ax.stem([x[0] for x in env], [x[1] for x in env], markerfmt = ' ')
        for x, y in env:
            if y > 0.01:
                ax.annotate('{}'.format(round(x,4)), xy = (x,y))

        plt.savefig('/Users/Aaron/local/envFinder/testFiles/test_env_{}.pdf'.format(i))



if __name__ == "__main__":
    main()
