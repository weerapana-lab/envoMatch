
from typing import List
import sys
import numpy as np

from .utils import find_nearest_index

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
                out.write(' -> {}, {}'.format(self.link.point.mz, self.link.point.int))
            out.write('\n')
        else:
            out.write('{}\t{}'.format(self.point.mz, self.point.int))
            if print_link:
                out.write('\t{}\t{}'.format(self.link.point.mz, self.link.point.int))
            out.write('\n')


def _getToleranceFunction(toleranceType: str, range: float):
    '''
    Get function to determine whether m/z matches are in the
    tolerance range.

    :param toleranceType: Tolerance units, either ppm or th.
    :type toleranceType: str
    :param range: Tolerance above and below.
    :type range: float
    :return: function
    '''

    if toleranceType == 'th':
        return lambda value, compare: abs(value - compare) <= range
    elif toleranceType == 'ppm':
        return lambda value, compare: abs(value - (compare * (range / 1e6))) <= value
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
                 tolerance: float = 50, toleranceType:str = 'ppm'):
        self._inRange = _getToleranceFunction(toleranceType, tolerance)
        self._actual : List = actual
        self._theoretical : List = theoretical
        self.initialized = False


    def setActual(self, actual: list):
        self._clearLinks(self._theoretical)
        self.initialized = False
        self._actual = actual


    def setTheoretical(self, theoretical: List):
        self._clearLinks(self._actual)
        self.initialized = False
        self._theoretical = theoretical


    def iterActual(self):
        return enumerate(self._actual)


    def iterTheoretical(self):
        return enumerate(self._theoretical)


    def annotate(self, remove_unlabeled: bool = False, normalize = True):

        self._clearLinks(self._theoretical)
        self._clearLinks(self._actual)

        for i, peak in enumerate(self._theoretical):
            best_i = find_nearest_index(self._actual, peak,
                                        valueKey= lambda x: x.point.mz,
                                        arrKey = lambda x: x.point.mz)

            if self._inRange(self._actual[best_i].point.mz, peak.point.mz):
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
        cor = np.corrcoef(x = [x.point.int for x in self._actual],
                          y = [x.link.point.int if x.link is not None else 0 for x in self._actual])
        self.envScore = cor[0,1]


    def plotEnv(self, fname:str):
        pass


