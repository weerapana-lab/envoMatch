
import re

_ATOM_RE = r'(\(\d+\))?([A-Z][a-z]?)(\d+)?'

def inRange(value, compare, range):
    '''
    Determine whether the difference between value and compare is <= range
    :param value: value
    :type numeric:
    :param compare: value to compare to
    :type numeric:
    :param range: ± range
    :type numeric:
    :return: True if value and compare are in range
    '''
    return abs(value - compare) <= range


def _numericalCompare(lhs, rhs,
                      lkey = lambda x: x,
                      rkey = lambda x: x):
    if lkey(lhs) < rkey(rhs):
        return -1
    elif lkey(lhs) > rkey(rhs):
        return 1
    else:
        return 0


def lower_bound(array, value,
                arrKey = lambda x: x,
                valueKey = lambda x: x,
                comp = _numericalCompare):
    '''
    Return the value i such that all e in array[:i] have e < value
    and all e in array[i:] have e >= value.

    :param array: sorted, iterable object
    :param value: numeric value (float or int)
    :param arrKey: key to access compare value of elements in array
    :param valueKey: key to access compare value in value
    :param comp: compare method
    :return: index of closest value in array
    :type int:
    '''

    lo = 0
    hi = len(array) - 1

    while lo < hi:
        mid = (lo+hi)//2

        comp_temp = comp(array[mid], value, arrKey, valueKey)
        if comp_temp == -1:
            lo = mid + 1
        elif comp_temp == 1 or comp_temp == 0:
            hi = mid
        else:
            raise RuntimeError('Invalid compare value!')

    return lo


def find_nearest_index(array, value,
                       arrKey = lambda x: x,
                       valueKey = lambda x: x,
                       comp = _numericalCompare):
    '''
    Get the index of the closest value to value in array.

    :param array: iterable object
    :param value: numeric value (float or int)
    :param arrKey: key to access compare value of elements in array
    :param valueKey: key to access compare value in value
    :param comp: compare method
    :return: index of closest value in array
    :type int:
    '''

    if len(array) == 0 or comp(value, array[0], valueKey, arrKey) == -1:
        return 0
    if comp(value, array[-1], valueKey, arrKey) == 1:
        return len(array) - 1

    lo  = 0
    hi = len(array) - 1

    while lo <= hi:
        mid: int = int((hi + lo) / 2)
        comp_temp: int = comp(value, array[mid], valueKey, arrKey)

        if comp_temp == -1:
            hi = mid - 1
        elif comp_temp == 1:
            lo = mid + 1
        elif comp_temp == 0:
            return mid
        else:
            raise RuntimeError('Invalid compare value!')

    return lo if (arrKey(array[lo]) - valueKey(value)) <\
                 (valueKey(value) - arrKey(array[hi])) else hi


def getRank(lst, lstKey = lambda x: x):
    '''
    Get a list of values in lst ranked by lstKey.

    :param lst: List of numeric values or objects
    :param lstKey: Accessor for value to order elements in lst by
    :return: A list or ranks of values in lst in the original order given
    '''
    values = [(i, lstKey(x)) for i, x in enumerate(lst)]
    values = sorted(values, key = lambda x: x[1], reverse = True)
    return [x[0] for x in values]


def formula_to_pyteomics_formula(formula: str):
    '''
    Convert a properly formatted formula to the stupid format used
    by the pyteomics package.

    :param formula: Properly formatted formula.
    :return: pyteomics format formula.
    '''
    ret = ''
    for atom in re.findall(_ATOM_RE, formula):
        ret += '{}{}{}'.format(atom[1],
                               '' if atom[0] == '' else '[{}]'.format(atom[0]),
                               atom[2])

    return ret



