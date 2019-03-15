
import utils

#def norm_intensity()


def findEnvelope(theorEnv, spec):

    #iterate through ions in theorEnv and find nearest match in sped
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
