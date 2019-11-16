
from typing import List
import sys
import numpy as np
import matplotlib.pyplot as plt

from .utils import lower_bound, inRange

class Isotope(object):
    __slots__ = ['mz', 'int']

    def __init__(self, mz: float = 0,
                 intensity: float = 0):
        self.mz = mz
        self.int = intensity


class DataPoint(object):

    def __init__(self, mz, intensity):
        self.point = Isotope(mz, intensity)
        self.link = None

    def found(self) -> bool:
        return self.link == None

    def print(self, out = sys.stdout,
              pretty_print = True, print_link = True):
        '''
        Print DataPoint to out.

        :param out: Stream to print to.
        :param pretty_print: If true pretty formatted text will be printed,
        if false, easily parsed text is printed.
        :type pretty_print: bool
        :param print_link: Should the mz and int of linked point be printed?
        :type print_link: bool
        '''

        if pretty_print:
            out.write('{}, {}'.format(self.point.mz, self.point.int))
            if print_link:
                if self.link is not None:
                    out.write(' -> {}, {}'.format(self.link.point.mz, self.link.point.int))
            out.write('\n')
        else:
            out.write('{}\t{}'.format(self.point.mz, self.point.int))
            if print_link:
                if self.link is not None:
                    out.write(' -> {}, {}'.format(self.link.point.mz, self.link.point.int))
            out.write('\n')


def _getTolerance(toleranceType: str, _range: float):
    '''
    Get function to return m/z +- tolerance.

    :param toleranceType: Tolerance units, either ppm or th.
    :type toleranceType: str
    :param _range: Tolerance above and below.
    :type _range: float
    :return: function
    '''

    if toleranceType == 'th':
        return lambda value : _range
    elif toleranceType == 'ppm':
        return lambda value: (value * (_range / 1e6))
    else:
        raise ValueError('{} is a invalid argument for toleranceType'.format(toleranceType))


class ConsensusEnvelope(object):

    @staticmethod
    def _clearLinks(lst: List):
        for i, _ in enumerate(lst):
            lst[i].link = None

    #@staticmethod
    #def _iterate(iterable):
    #    for i in

    def __init__(self, actual: List[DataPoint] = None, theoretical: List[Isotope] = None,
                 tolerance: float = 50, toleranceType:str = 'ppm',
                 best_match_tie:str = 'intensity', envScoreCuttoff:float = 0.8,
                 sequence: str = None):
        self.getTolerance = _getTolerance(toleranceType, tolerance)
        self._actual : List = actual
        self._theoretical : List = theoretical
        self._mono_mz = None
        self._mono_index = None
        self._best_match_tie = best_match_tie
        self.initialized = False
        self.envScoreCuttoff = envScoreCuttoff
        self.sequence = sequence


    def setActual(self, actual: list):
        self._clearLinks(self._theoretical)
        self.initialized = False
        self._actual = actual


    def setTheoretical(self, theoretical: List):
        self._clearLinks(self._actual)
        self.initialized = False
        self._theoretical = theoretical


    def set_mono(self, mono_mz: float = None,
                 mono_mass: float = None, charge: int = None):

        if mono_mz is not None:
            self._mono_mz = mono_mz
        elif mono_mass is not None and charge is not None:
            self._mono_mz = (mono_mass + charge) / charge
        else:
            raise RuntimeError('Either mono_mz or mono_mass and charge are required!')

        #indices = [x for x in self._theoretical if self._inRange(x.point.mz, self._mono_mz)]
        indices = list()
        for i, x in enumerate(self._theoretical):
            if inRange(x.point.mz, self._mono_mz, self.getTolerance(self._mono_mz)):
                indices.append(i)

        if len(indices) != 1:
            raise RuntimeError('Could not find mono_mz in theoretical env!')
        self._mono_index = indices[0]


    def iterActual(self):
        return enumerate(self._actual)


    def iterTheoretical(self):
        return enumerate(self._theoretical)


    def annotate(self, remove_unlabeled: bool = False, normalize = True):
        '''
        Annotate actual spectra with theoretical envelope.

        :param remove_unlabeled: Should unlabeled points in actual spectra be removed?
        :type remove_unlabeled: bool
        :param normalize: Should both theoretical and actual intensities be normalized to 1?
        :type normalize: bool
        :return:
        '''

        self._clearLinks(self._theoretical)
        self._clearLinks(self._actual)

        for i, peak in enumerate(self._theoretical):

            low_i = lower_bound(self._actual, peak,
                                valueKey = lambda x: x.point.mz - self.getTolerance(x.point.mz),
                                arrKey = lambda x: x.point.mz)

            range_list = list()
            #best_i = 0
            for j in range(low_i, len(self._actual), 1):
                if self._actual[j].point.mz > peak.point.mz + self.getTolerance(peak.point.mz):
                    break
                if inRange(self._actual[j].point.mz, peak.point.mz, self.getTolerance(self._actual[j].point.mz)):
                    range_list.append(j)

            if len(range_list) == 0:
                continue
            elif len(range_list) == 1:
                best_i = range_list[0]
            else:
                tempMax = range_list[0]
                if self._best_match_tie == 'intensity':
                    def _temp_compare(lhs, rhs):
                        return lhs.point.int > rhs.point.int
                else:
                    def _temp_compare(lhs, rhs):
                        return lhs.point.mz > rhs.point.mz

                for j in range_list:
                    if _temp_compare(self._actual[j], self._actual[tempMax]):
                        tempMax = j
                best_i = tempMax

            self._actual[best_i].link = peak
            self._theoretical[i].link = self._actual[best_i]

        if remove_unlabeled:
            self._actual = [x for x in self._actual if x.link is not None]

        if normalize:
            try:
                max_int = max([x.point.int for x in self._actual if x.link is not None])
            except ValueError:
                sys.stdout.write('No points in envelope found!\n')
                max_int = max(self._actual, key = lambda x: x.point.int).point.int

            for i in range(len(self._actual)):
                self._actual[i].point.int /= max_int

            max_int = max(self._theoretical, key = lambda x: x.point.int).point.int
            for i in range(len(self._theoretical)):
                self._theoretical[i].point.int /= max_int

        self.calcEnvScore()
        self.initialized = True


    def calcEnvScore(self):
        cor = np.corrcoef(x = [x.point.int for x in self._theoretical],
                          y = [x.link.point.int if x.link is not None else 0 for x in self._theoretical])
        self.envScore = cor[0,1]


    def plotEnv(self, ax, title: str = None, isBest:bool = False):

        #general properties of plot
        ax.margins(y = 0)
        ax.set_ylim(0, 1.1)
        ax.set_ylabel('Relative intensity')

        _title = None
        if title is None:
            if self.sequence is not None:
                _title = '{}, r = {:.2f}'.format(self.sequence, self.envScore)
        else:
            _title = title
        if _title is not None:
            ax.set_title(_title, {'fontweight': 'bold' if isBest else plt.rcParams['axes.titleweight']}, 'center')

        #plot actual spectra
        _, stemlines, _ = ax.stem([x.point.mz for x in self._actual],
                                   [x.point.int for x in self._actual],
                                   basefmt = ' ', markerfmt = ' ')

        for i in range(len(stemlines)):
            if self._actual[i].link is None:
                plt.setp(stemlines[i], color = 'black', alpha = 0.4)
                #plt.setp(markerlines[i], 'Visible', False)
            else:
                plt.setp(stemlines[i], 'color', 'blue')

        #plot theoretical envelope
        for i, it in enumerate(self._theoretical):
            fill = 'white'
            if self._mono_index is not None:
                if i == self._mono_index:
                    fill = 'green' if it.link is not None else 'red'

            outline = 'green' if it.link is not None else 'red'

            markerlines, _, _ = ax.stem([it.point.mz], [0], basefmt=' ', linefmt=' ', markerfmt='D')
            plt.setp(markerlines, mfc = fill, mec = outline)


        ax.plot([x.point.mz for x in self._theoretical],
                [x.point.int for x in self._theoretical],
                'o--{}'.format('g' if self.envScore > self.envScoreCuttoff else 'r'),
                alpha = 0.6)

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)







