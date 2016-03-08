from builtins import zip
from builtins import range
def enum(*sequential, **named):
    """creates an Enum type with given values"""
    enums = dict(list(zip(sequential, list(range(len(sequential))))), **named)
    enums['reverse'] = dict((value, key) for key, value in list(enums.items()))
    return type('Enum', (object, ), enums)