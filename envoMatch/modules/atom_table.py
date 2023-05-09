
import sys
from typing import Dict
from collections import Counter
from sortedcontainers import SortedList

from pyteomics.mass import Composition
from pyteomics.mass.mass import isotopologues

from .utils import inRange, find_nearest_index


DEFAULT_COMPOSITIONS = {'A': Counter({'C': 3, 'H': 5, 'O': 1, 'N': 1}),
                        'C': Counter({'C': 5, 'H': 8, 'O': 2, 'N': 2, 'S': 1}),
                        'D': Counter({'C': 4, 'H': 5, 'O': 3, 'N': 1}),
                        'E': Counter({'C': 5, 'H': 7, 'O': 3, 'N': 1}),
                        'F': Counter({'C': 9, 'H': 9, 'O': 1, 'N': 1}),
                        'G': Counter({'C': 2, 'H': 3, 'O': 1, 'N': 1}),
                        'H': Counter({'C': 6, 'H': 7, 'O': 1, 'N': 3}),
                        'I': Counter({'C': 6, 'H': 11, 'O': 1, 'N': 1}),
                        'K': Counter({'C': 6, 'H': 12, 'O': 1, 'N': 2}),
                        'L': Counter({'C': 6, 'H': 11, 'O': 1, 'N': 1}),
                        'M': Counter({'C': 5, 'H': 9, 'O': 1, 'N': 1, 'S': 1}),
                        'N': Counter({'C': 4, 'H': 6, 'O': 2, 'N': 2}),
                        'P': Counter({'C': 5, 'H': 7, 'O': 1, 'N': 1}),
                        'Q': Counter({'C': 5, 'H': 8, 'O': 2, 'N': 2}),
                        'R': Counter({'C': 6, 'H': 12, 'O': 1, 'N': 4}),
                        'S': Counter({'C': 3, 'H': 5, 'O': 2, 'N': 1}),
                        'T': Counter({'C': 4, 'H': 7, 'O': 2, 'N': 1}),
                        'V': Counter({'C': 5, 'H': 9, 'O': 1, 'N': 1}),
                        'W': Counter({'C': 11, 'H': 10, 'O': 1, 'N': 2}),
                        'Y': Counter({'C': 9, 'H': 9, 'O': 2, 'N': 1}),
                        'U': Counter({'C': 5, 'H': 8, 'O': 2, 'N': 2, 'Se': 1}),
                        'C_term': Counter({'H': 1, 'O': 1}),
                        'N_term': Counter({'H': 1}),
                        '*': Counter({'H': -1, 'O': 1, 'N': -1})}


class AtomTable:
    _nTermStr = 'N_term'
    _cTermStr = 'C_term'

    _atom_masses = {"C": 12,
                    "H": 1.00783,
                    "H+": 1.00783,
                    "O": 15.99491,
                    "N": 14.00307,
                    "S": 31.97207,
                    "P": 30.97376,
                    "N[15]": 15.00011,
                    "H[2]": 2.0141,
                    "C[13]": 13.00335,
                    "Se": 79.91652,
                    "Cl": 34.96885,
                    "Br": 78.91834}

    def __init__(self, fname: str=None):
        self.fname = fname
        self.compositions = dict()

    def _cleanLine(self, line):
        elems = [y for y in [x.strip() for x in line.split('\t')] if y and y[:1] != ';']
        return elems[0], elems[1:]


    def _read(self):
        if self.fname is None:
            self.compositions = DEFAULT_COMPOSITIONS
            return

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
                self.compositions[elems[0]] = Counter({k:int(v) for k, v in zip(headers[1:], elems[1:]) if v != '0'})
                continue


    def read(self, fname: str=None):

        #check fname
        if fname is not None:
            self.fname = fname
        try:
            self._read()
        except (IOError, RuntimeError) as e:
            sys.stderr.write('{}\n'.format(e))
            return False
        except Exception as e:
            sys.stderr.write('An unknown error occurred.\n{}\n'.format(e))
            return False

        return True


    def _getComposition(self, seq: str, charge=0,
                       nTerm=True, cTerm=True):

        tempDict = Counter()
        for aa in seq:
            tempDict.update(self.compositions[aa])

        if nTerm:
            tempDict.update(self.compositions[AtomTable._nTermStr])
        if cTerm:
            tempDict.update(self.compositions[AtomTable._cTermStr])

        if charge != 0:
            tempDict['H+'] = charge

        return tempDict


    def getComposition(self, seq: str, charge: int=1,
                       nTerm=True, cTerm=True):

        return Composition(self._getComposition(seq, charge=charge, nTerm=nTerm, cTerm=cTerm))


    def getMass(self, seq: str=None, composition: Dict=None,
                charge: int=0, nTerm=True, cTerm=True) -> float:
        '''
        Calculate mass of sequence or formula.

        Function accepts `seq` xor `composition`.

        :param seq: peptide sequence
        :param composition: Peptide composition.
        :param charge: Peptide charge.
        If None, 'H+' in composition is used, else specified charge is used.
        '''

        if sum([1 for x in [seq, composition] if x is not None]) != 1:
            raise RuntimeError('Either seq xor composition must be specified')

        if seq is not None:
            _comp = self.getComposition(seq, nTerm, cTerm, 0 if charge is None else charge)
        else:
            _comp = composition
            if charge is not None:
                _comp['H+'] = charge

        ret = float(0)
        for atom, count in _comp.items():
            ret += self._atom_masses[atom] * count

        return ret

class AverageIsotope():
    __slots__ = ['_count', '_mzSum', 'intSum', 'avg_mz']

    def __init__(self, mzSum: float=0, intSum: float=0):
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
                threshold: float=1e-3,
                mass_defect_match_range: float=0.05,
                combineDefects=True) -> list:

    #get all isotopologues with reasonable abundance
    temp = list()
    for (c, a) in isotopologues(composition,
                                report_abundance=True,
                                isotope_threshold=5e-3,
                                overall_threshold=threshold):
        charge = c.pop('H+[1]')
        temp.append((c.mass(charge=charge, charge_carrier='H+'), a))

    #combine and average mass defects
    if combineDefects:
        masses = SortedList(key=lambda x: x.avg_mz)
        for isotope in temp:
            idx = find_nearest_index(masses, isotope[0], arrKey=lambda x: x.avg_mz)

            if masses and inRange(isotope[0], masses[idx].avg_mz, mass_defect_match_range):
                masses[idx].append(isotope)
            else:
                masses.add(AverageIsotope(isotope[0], isotope[1]))

        return [x.value() for x in masses]
    return temp

