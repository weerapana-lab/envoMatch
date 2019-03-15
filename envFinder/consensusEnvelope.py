
import typing

import utils

class Isotope(object):
    __slots__ = ['mz', 'intensity', 'rank' 'found']

    def __init(self, mz: float = 0,
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
        self.name = name
        self.isotopes: ConsensusIsotope.IsotopesType = {name: isotope}

    #def addSample(self, sample):

    #adding and inserting functions


def _getToleranceFunction(toleranceType: str):
    if toleranceType == 'th':
        return lambda value, compare, range: abs(value - compare) <= range
    elif toleranceType == 'ppm':
        return lambda value, compare, range: abs(value - (compare * (range / 1e6))) <= value
    else:
        raise ValueError('{} is a invalid argument for toleranceType'.format(toleranceType))


class ConsensusEnvelope(object):

    def __init__(self,
                 isotopeTolerance: float = 0.02,
                 toleranceType: str = 'th'):

        self.envelope: typing.List[ConsensusIsotope] = list()
        self.isotopeTolerance: float = isotopeTolerance
        self.inRange = _getToleranceFunction(toleranceType)

        #peptide data, maybe inhereted from a Peptide base class


    def initalize(self, name: str,
                  envelope: typing.List[typing.Tuple[float, float]]):

        ranks = utils.getRank(envelope, lstKey = lambda x: x[1])
        for i, isotope in enumerate(envelope):
            self.envelope.append(ConsensusIsotope(isotope[0], Isotope(mz = isotope[0],
                                                                      intensity = isotope[1],
                                                                      rank = ranks[i],
                                                                      found = True)))


    #def getEnvelope(sampleName)
