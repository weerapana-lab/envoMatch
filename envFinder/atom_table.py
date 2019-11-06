
import sys
from collections import Counter
from sortedcontainers import SortedList

from pyteomics.mass import Composition
from pyteomics.mass.mass import isotopologues

import utils

class AtomTable:
    _nTermStr = 'N_term'
    _cTermStr = 'C_term'

    def __init__(self, fname:str = ''):
        self.fname = fname
        self.composition = Counter()

    def _cleanLine(self, line):
        elems = [y for y in [x.strip() for x in line.split('\t')] if y and y[:1] != ';']
        return elems[0], elems[1:]


    def _read(self):
        with open(self.fname, 'rU') as inF:
            lines = inF.readlines()

        headers = dict()
        headersLen = 0
        for line in lines:
            tag, elems = self._cleanLine(line)
            if elems[0][0] == '#':
                continue

            if tag == 'H':
                headers = [x for x in elems]
                headersLen = len(headers)
                continue
            elif tag == 'R':
                if not headers or len(elems) != headersLen:
                    raise RuntimeError('Badly formed atom_table file!')
                self.composition[elems[0]] = Counter({k:int(v) for k, v in zip(headers[1:], elems[1:]) if v != '0'})
                continue


    def read(self, fname: str = ''):

        #check fname
        if not self.fname:
            self.fname = fname
        if not self.fname:
            sys.stderr.write('fname is required!\n')
            return False

        try:
            self._read()
        except (IOError, RuntimeError) as e:
            sys.stderr.write('{}\n'.format(e))
            return False
        except Exception as e:
            sys.stderr.write('An unknown error occured.\n{}\n'.format(e))
            return False

        return True


    def getComposition(self, seq: str, charge: int = 1,
                       nTerm = True, cTerm = True):

        tempDict = Counter()
        for aa in seq:
            tempDict.update(self.composition[aa])

        if nTerm: tempDict.update(self.composition[AtomTable._nTermStr])
        if cTerm: tempDict.update(self.composition[AtomTable._cTermStr])

        return Composition(tempDict, charge = charge)


class AverageIsotope():
    __slots__ = ['_count', '_mzSum', 'intSum', 'avg_mz']

    def __init__(self, mzSum: float = 0, intSum: float = 0):
        self._count = int(0) if mzSum == 0 else 1
        self._mzSum = mzSum
        self.intSum = intSum
        self.avg_mz = mzSum

    def append(self, isotope: tuple):
        '''
        Append isotope

        :param isotope: (mz, intensity)
        :type tuple:
        '''
        self._count += 1
        self._mzSum += isotope[0]
        self.intSum += isotope[1]
        self.avg_mz = self._mzSum / self._count

    def value(self) -> tuple:
        '''
        Get tuple of (average mz, intensity sum)
        :return:
        '''
        return (self._mzSum / self._count, self.intSum)


    def __lt__(self, rhs: float) -> bool:
        return self.avg_mz < rhs

    def __gt__(self, rhs: float) -> bool:
        return self.avg_mz > rhs

    def compare(self, rhs) -> int:
        if self.avg_mz < rhs.avg_mz:
            return -1
        elif self.avg_mz > rhs.avg_mz:
            return 1
        else:
            return 0


def getEnvelope(composition: Composition,
                threshold: float = 1e-3,
                mass_defect_match_range: float = 0.05,
                combineDefects = True) -> list:

    #get all isotopologues with reasonable abundance
    temp = list()
    for (c, a) in isotopologues(composition,
                                report_abundance=True,
                                isotope_threshold=5e-3,
                                overall_threshold=threshold):
        c['H+'] = c.pop('H+[1]')
        temp.append((c.mass(), a))

    #combine and average mass defects
    if combineDefects:
        masses = SortedList(key = lambda x: x.avg_mz)
        for isotope in temp:
            idx = utils.find_nearest_index(masses, isotope[0], arrKey = lambda x: x.avg_mz)

            if masses and utils.inRange(isotope[0], masses[idx].avg_mz, mass_defect_match_range):
                masses[idx].append(isotope)
            else:
                masses.add(AverageIsotope(isotope[0], isotope[1]))

        return [x.value() for x in masses]
    else:
        return temp
