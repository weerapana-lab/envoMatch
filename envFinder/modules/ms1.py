
import sys
from pyteomics import ms1, mzml, mzxml

_MZ_KEY = 'mz'
_INT_KEY = 'int'

class Ms1File(object):

    @staticmethod
    def getMzRange(mono_mz, charge, range = 10):
       return ((mono_mz - range) / charge, (mono_mz + range) / charge)

    @staticmethod
    def getReadFxn(file_type):
        if file_type == 'ms1':
            return ms1.read
        elif file_type == 'mzXML':
            return mzxml.read
        elif file_type == 'mzML':
            return mzml.read
        else:
            raise RuntimeError('{} is an invalid file type!'.format(file_type))

    def _getIDStr(self, id):
        if self.file_type == 'ms1':
            return str(id).zfill(6)
        else:
            return str(id)

    def __init__(self, fname = None, file_type = 'mzXML'):
        self.file_type = file_type
        if fname is None:
            self.fname = str()
        else: self.read(fname, file_type)

    def read(self, fname, file_type = 'mzXML'):
        self.fname = fname
        _read = Ms1File.getReadFxn(file_type)
        self.dat = _read(self.fname, use_index = True)

    def getSpectra(self, scan, mz_range):
        '''
        Return spectra at scan in the mz_range

        :param scan: Scan number to fetch
        :type scan: int
        :param mz_range: mz range to return. If None, the entire scan is returned.
        :type mz_range: Tuple(float, float)
        :return: Dict with arrays for mz and int
        '''

        vals = {'m/z array': _MZ_KEY, 'intensity array': _INT_KEY}
        try:
            spec = self.dat.get_by_id(self._getIDStr(scan))
        except KeyError as e:
            sys.stderr.write('Scan ID: {} not found!\n'.format(scan))

        if mz_range is None:
            selection = [True for _ in spec['m/z array']]
        else:
            selection = list(map(lambda x, y: x and y,
                                 spec['m/z array'] >= mz_range[0],
                                 spec['m/z array'] <= mz_range[1]))

        return {v: spec[k][selection] for k, v in vals.items()}









