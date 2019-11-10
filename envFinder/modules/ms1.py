
import pdb
import numpy as np
from typing import Dict, Tuple, List

from pyteomics import ms1

#from .utils import find_nearest_index

_MZ_KEY = 'mz'
_INT_KEY = 'int'

class Ms1File(object):

    @staticmethod
    def getMzRange(mono_mz: float, charge: int, range : int = 10) -> Tuple:
       return ((mono_mz - range) / charge, (mono_mz + range) / charge)

    def __init__(self, fname: str = None):
        if fname is None:
            self.fname = str()
        else: self.read(fname)

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

        return {v: spec[k][selection] for k, v in vals.items()}









