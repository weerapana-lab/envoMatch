
import sys
import re
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

    def _get_scan(self, scan_id):
        if self.file_type == 'ms1':
            return self.dat.get_by_id(str(scan_id).zfill(6))
        elif self.file_type == 'mzML':
            return self.dat[self._scan_map[scan_id]]
        return self.dat.get_by_id(str(scan_id))

    def __init__(self, fname=None, file_type='mzXML', **kwargs):
        '''
        Default constructor.

        Parameters
        ----------
        fname: str
            Path to file to read. If None, no file is read.
        file_type: str
            Input file type. One of (mzXML, mzML, ms1)
        **kwargs
            Additional arguments passes to self.read
        '''

        self.file_type = file_type
        self.precursors = None
        self._scan_map = None
        if fname is None:
            self.fname = str()
        else: self.read(fname, file_type, **kwargs)

    def _build_precursor_list(self):
        if self.file_type == 'mzXML':
            _scans = {int(x['num']):x['msLevel'] for x in self.dat}
        elif self.file_type == 'mzML':
            scanPattern = re.compile(r'scan=([0-9]+)')

            # test scan pattern on the first scan
            match = scanPattern.search(self.dat[0]['id'])
            if match:
                use_index = False
            else:
                use_index = True
                sys.stderr.write('WARN: Failed to match scan for file {}.\n\tUsing index instead.\n')

            self._scan_map = dict()
            _scans = dict()
            for i, s in enumerate(self.dat):
                if use_index:
                    _scans[s['index'] + 1] = s['ms level']
                    self._scan_map[s['index']] = i
                else:
                    match = scanPattern.search(s['id'])
                    if match:
                        _scan = int(match.group(1))
                        _scans[_scan] = s['ms level']
                        self._scan_map[_scan] = i
                    else:
                        raise RuntimeError('Failed to find scan number for line {} in file {}'.format(s['id'], self.fname))
        else:
            assert(False)

        self.precursors = dict()
        for scan in sorted(_scans.keys()):
            if _scans[scan] == 1:
                pre_scan = scan
            elif _scans[scan] == 2:
                self.precursors[scan] = pre_scan

    def read(self, fname, file_type='mzXML', build_precursor_list=False):
        '''
        Read ms1 file.

        Parameters
        ----------
        fname: str
            Name of file to read.
        file_type: str
            File type. One of (mzXML, mzML, ms1)
        build_precursor_list: bool
            Should self.precursors be populated?
        '''

        self.fname = fname
        _read = Ms1File.getReadFxn(file_type)
        self.dat = _read(self.fname, use_index=True)
        if build_precursor_list:
            self._build_precursor_list()

    def get_precursor_scan(self, scan):
        '''
        Get precursor scan for `scan`.

        Raises
        ------
        KeyError
            If `scan` does not exist in ms file.
        '''
        try:
            return self.precursors[scan]
        except KeyError as e:
            raise KeyError('Scan: {} does not exist in ms0 file: {}'.format(scan, self.fname))

    def get_spectra(self, scan, mz_range):
        '''
        Return spectra at scan in the mz_range

        Parameters
        ----------
        scan: int
            Scan number to fetch
        mz_range: Tuple(float, float)
            mz range to return. If None, the entire scan is returned.

        Returns
        -------
            Dict with arrays for mz and int
        '''

        vals = {'m/z array': _MZ_KEY, 'intensity array': _INT_KEY}
        try:
            spec = self._get_scan(scan)
        except KeyError as e:
            sys.stderr.write('Scan ID: {} not found!\n'.format(scan))
            return None

        if mz_range is None:
            selection = [True for _ in spec['m/z array']]
        else:
            selection = list(map(lambda x, y: x and y,
                                 spec['m/z array'] >= mz_range[0],
                                 spec['m/z array'] <= mz_range[1]))

        return {v: spec[k][selection] for k, v in vals.items()}

