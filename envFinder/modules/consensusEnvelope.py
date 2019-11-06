
from typing import List
from sortedcontainers import SortedList

import utils

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

    def __init__(self, actual: list = None, theoretical = None,
                 tolerance: float = 50, toleranceType:str = 'ppm'):
        self._inRange = _getToleranceFunction(toleranceType, tolerance)
        self._actual : List = actual
        self._theoretical : List = theoretical
        #self._all_labels : List = list()
        self.initialized = False


    def setActual(self, actual: list):
        self._clearLinks(self._theoretical)
        #self._all_labels = None
        self.initialized = False
        self._actual = actual

    def setTheoretical(self, theoretical: List):
        self._clearLinks(self._actual)
        #self._all_labels = None
        self.initialized = False
        self._theoretical = theoretical

    def iterActual(self):
        return enumerate(self._actual)

    def iterTheoretical(self):
        return enumerate(self._theoretical)

    #def iterLabels(self):
    #    return enumerate(self._all_labels)

    def annotate(self):

        self._clearLinks(self._theoretical)
        self._clearLinks(self._actual)

        for i, peak in enumerate(self._theoretical):
            best_i = utils.find_nearest_index(self._actual, peak,
                                         arrKey = lambda x: x.point.mz)

            if self._inRange(self._actual[best_i].point.mz, peak):
                self._actual[best_i].link = peak
                self._theoretical[i].link = self._actual[best_i].link

        self.initialized = True

    def calcEnvScore(self):
        pass

    def plotEnv(self, fname:str):
        pass


