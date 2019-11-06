
import typing
from sortedcontainers import SortedList

import utils

class Isotope(object):
    __slots__ = ['mz', 'intensity', 'rank', 'found']

    def __init__(self, mz: float = 0,
               intensity: float = 0,
               rank: int = 0,
               found: bool = False):
        self.mz = mz
        self.intensity = intensity
        self.rank = rank
        self.found = found


class ConsensusIsotope(object):
    IsotopesType = typing.NewType('IsotopesType', typing.Dict[str, Isotope])

    def __init__(self, mz: float, isotope: Isotope, name: str = 'theoretical'):
        self.mz =  mz
        self.mz_sum = mz
        self.count = int(0) if mz == 0 else 1
        self.avg_mz = mz
        self._isotopes: ConsensusIsotope.IsotopesType = {name: isotope}


    def addSample(self, mz: float, isotope: Isotope, name: str = 'experimental'):
        self.mz_sum += mz
        self.count += 1
        self.avg_mz = self.mz_sum / self.count

        self._isotopes[name] = isotope


    #adding and inserting functions


def _getToleranceFunction(toleranceType: str, range: float):
    if toleranceType == 'th':
        return lambda value, compare: abs(value - compare) <= range
    elif toleranceType == 'ppm':
        return lambda value, compare: abs(value - (compare * (range / 1e6))) <= value
    else:
        raise ValueError('{} is a invalid argument for toleranceType'.format(toleranceType))


class ConsensusEnvelope(object):

    DefaultTheoreticalName = 'theoretical'
    DefaultExperimentalName = 'experimental'

    def __init__(self,
                 isotopeTolerance: float = 0.02,
                 toleranceType: str = 'th'):

        self.envelope = SortedList(key = lambda x: x.mz)
        #self.isotopeTolerance: float = isotopeTolerance
        self.inRange = _getToleranceFunction(toleranceType, isotopeTolerance)
        self.experimentalNames = list()
        self.theoreticalName = ConsensusEnvelope.DefaultTheoreticalName

        #peptide data, maybe inhereted from a Peptide base class


    def initalize(self,
                  envelope: typing.List[typing.Tuple[float, float]],
                  name: str = DefaultTheoreticalName):

        self.theoreticalName = name

        ranks = utils.getRank(envelope, lstKey = lambda x: x[1])
        for i, isotope in enumerate(envelope):
            self.envelope.add(ConsensusIsotope(isotope[0],
                                               Isotope(mz = isotope[0],
                                                       intensity = isotope[1],
                                                       rank = ranks[i],
                                                       found = True),
                                               name = name))


    def _addSample(self, name: str):
        for isotope in self.envelope:
            isotope._isotopes[name] = Isotope(found = False)


    def add(self,
            envelope: typing.List[typing.Tuple[float, float]],
            name: str = DefaultExperimentalName):

        self.experimentalNames.append(name)

        ranks = utils.getRank(envelope, lstKey=lambda x: x[1])
        for i, isotope in enumerate(envelope):
            idx = utils.find_nearest_index(self.envelope, isotope[0],
                                           arrKey = lambda x: x.mz)

            if self.inRange(self.envelope[idx].mz, isotope[0]):
                self.envelope[idx].addSample(isotope[0],
                                             Isotope(mz = isotope[0],
                                                     intensity = isotope[1],
                                                     rank = ranks[i],
                                                     found = True),
                                             name = name)
            else:
                self.envelope.add(ConsensusIsotope(isotope[0],
                                                   Isotope(found=False),
                                                   name = self.theoreticalName))

                for x in self.experimentalNames:
                    self._addSample(x)

                self.envelope.add(ConsensusIsotope(isotope[0],
                                                   Isotope(mz = isotope[0],
                                                           intensity=isotope[1],
                                                           rank = ranks[1],
                                                           found=False),
                                                   name=name))


    # def getEnvelope(self, ):


    #def getEnvelope(sampleName)
