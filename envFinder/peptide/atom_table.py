
import sys
from collections import Counter

from pyteomics.mass import Composition
from pyteomics.mass.mass import isotopologues

from .constatns import _nTermStr, _cTermStr

class AtomTable:
    def __init__(self, fname = ''):
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


    def read(self, fname = ''):

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


    def getComposition(self, seq, charge = 1,
                       nTerm = True, cTerm = True):

        tempDict = Counter()
        for aa in seq:
            tempDict.update(self.composition[aa])

        if nTerm: tempDict.update(self.composition[_nTermStr])
        if cTerm: tempDict.update(self.composition[_cTermStr])

        return Composition(tempDict, charge = charge)


def getEnvelope(composition, threshold = 1e-3):

    ret = list()
    for (c, a) in isotopologues(composition, report_abundance=True, isotope_threshold=5e-3,
                                overall_threshold=threshold):
        c['H+'] = c.pop('H+[1]')
        ret.append((c.mass(), a))

    return ret
