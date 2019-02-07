
import pyteomics
import sys
import operator
import matplotlib.pyplot as plt

import peptide
import parse_args

def main():
    args = parse_args.parseArgs()

    atomTableFname = '/Users/Aaron/local/envFinder/db/atom_tables/cit_diff_mod_atoms.txt'
    atomTable = peptide.AtomTable(atomTableFname)

    if not atomTable.read():
        exit()

    #for aa, comp in atomTable.composition.items():
    #    print('{} -> {}'.format(aa, comp))

    c = atomTable.getComposition('GERALDGERALDGERALDGERALD', 2)
    #atomTable.getComposition('GER*ALD')
    env = peptide.getEnvelope(c)

    for x in map(operator.itemgetter(0), env):
        print(x)

    fig, ax = plt.subplots()
    ax.stem([x[0] for x in env], [x[1] for x in env], markerfmt = ' ')
    #plt.show()
    plt.savefig('/Users/Aaron/local/envFinder/testFiles/test_env.pdf')

if __name__ == "__main__":
    main()
