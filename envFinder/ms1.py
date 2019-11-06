
import numpy as np
from typing import Dict, Tuple, List

from pyteomics import ms1

#from .utils import find_nearest_index

_MZ_KEY = 'mz'
_INT_KEY = 'int'

def findEnvelope(theorEnv, spec):

    #iterate through ions in theorEnv and find nearest match in spec
    temp = list()
    for mass in theorEnv:
        temp.append(utils.find_nearest_index(spec['m/z array'], mass[0]))

    temp = [(spec['m/z array'][x], spec['intensity array'][x]) for x in set(temp)]

    #sum intensity
    intSum = 0
    for ion in temp:
        intSum += ion[1]

    #convert intensities to a fraction of 1

    ret = list()
    for ion in temp:
        ret.append((ion[0], ion[1] / intSum))

    return ret


def annotateSpectra(spec: Dict, env: List, tolerance: float = 50, toleranceType:str = 'ppm') -> Dict:

   spec['found'] = np.array([False for _ in spec[_MZ_KEY]])

   inRange = _getToleranceFunction(toleranceType, tolerance)
   for peak in env:
       i = utils.find_nearest_index(spec[_MZ_KEY], peak)
       if inRange(spec[_MZ_KEY][i], peak):
           spec


class Ms1File(object):

    @staticmethod
    def getMzRange(mono_mz: float, charge: int, range : int = 10) -> Tuple:
       return ((mono_mz - range) / charge, (mono_mz + range) / charge)

    def __init__(self):
        self.fname = str()
        #self.dat = ms1()

    def read(self, fname: str):
        self.fname = fname
        self.dat = ms1.read(self.fname, use_index = True)

    def getSpectra(self, scan: int, mz_range : Tuple) -> Dict:
        '''
        Return spectra at scan in the mz_range

        :param scan: Scan number to fetch
        :type scan: int
        :param mz_range: mz range to return. If None, the entire scan is returned.
        :type mz_range: Tuple(float, float)
        :return: Dict with arrays for mz and int
        '''
        vals = {'m/z array': _MZ_KEY, 'intensity array': _INT_KEY}
        spec = self.dat.get_by_id(str(scan).zfill(6))
        if mz_range is None:
            selection = [True for _ in spec['m/z array']]
        else:
            selection = list(map(lambda x, y: x and y,
                                 spec['m/z array'] >= mz_range[0],
                                 spec['m/z array'] <= mz_range[1]))

        return {v: spec[selection] for k, v in vals.items()}









