def enum(*sequential, **named):
    """creates an Enum type with given values"""
    enums = dict(zip(sequential, range(len(sequential))), **named)
    enums['reverse'] = dict((value, key) for key, value in enums.items())
    return type('Enum', (object, ), enums)