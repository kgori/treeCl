"""
Data structures to hold model parameters and attributes such as alignment file paths
"""
from __future__ import absolute_import
from builtins import object
import json
import sys
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def setter_helper(fn, value):
    """
    Some property setters need to call a function on their value -- e.g. str() or float() -- before setting
    the property. But, they also need to be able to accept None as a value, which throws an error.
    This helper solves this problem

    :param fn: function to apply to value
    :param value: value to pass to property setter
    :return: result of fn(value), or None
    """
    return value if value is None else fn(value)

class BaseParameters(object):
    __slots__ = []

    @property
    def dict(self):
        return dict((item.lstrip('_'), getattr(self, item)) for item in self.__slots__)

    def read(self, fileobj):
        return json.load(fileobj)

    def construct_from_dict(self, dict):
        for k, v in dict.items():
            if '_{}'.format(k) in self.__slots__:
                setattr(self, k, v)

    def write(self, fileobj=sys.stdout, indent=None):
        json.dump(self.dict, fileobj, indent=indent)


class PartitionParameters(BaseParameters):
    __slots__ = ['_alpha', '_distances', '_frequencies', '_model', '_name', '_rates', '_variances']

    def __init__(self):
        self._alpha = None
        self._distances = None
        self._frequencies = None
        self._model = None
        self._name = None
        self._rates = None
        self._variances = None

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @property
    def distances(self):
        return self._distances

    @distances.setter
    def distances(self, value):
        self._distances = value

    @property
    def frequencies(self):
        return self._frequencies

    @frequencies.setter
    def frequencies(self, value):
        self._frequencies = value

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        self._model = setter_helper(str, value)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    @property
    def rates(self):
        return self._rates

    @rates.setter
    def rates(self, value):
        self._rates = value

    @property
    def variances(self):
        return self._variances

    @variances.setter
    def variances(self, value):
        self._variances = value


class Parameters(BaseParameters):
    __slots__ = ['_filename', '_likelihood', '_ml_tree', '_ms_tree', '_nj_tree', '_partitions', '_sse']

    def __init__(self):
        self._filename = None
        self._likelihood = None
        self._ml_tree = None
        self._ms_tree = None
        self._nj_tree = None
        self._partitions = []
        self._sse = None

    @property
    def filename(self):
        return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = setter_helper(str, value)

    @property
    def likelihood(self):
        return self._likelihood

    @likelihood.setter
    def likelihood(self, value):
        self._likelihood = setter_helper(float, value)

    @property
    def sse(self):
        return self._sse

    @sse.setter
    def sse(self, value):
        self._sse = setter_helper(float, value)

    @property
    def ml_tree(self):
        return self._ml_tree

    @ml_tree.setter
    def ml_tree(self, value):
        self._ml_tree = setter_helper(str, value)

    @property
    def nj_tree(self):
        return self._nj_tree

    @nj_tree.setter
    def nj_tree(self, value):
        self._nj_tree = setter_helper(str, value)

    @property
    def ms_tree(self):
        return self._ms_tree

    @ms_tree.setter
    def ms_tree(self, value):
        self._ms_tree = setter_helper(str, value)

    @property
    def partitions(self):
        if len(self._partitions) == 1:
            return self._partitions[0]
        else:
            return self._partitions

    @partitions.setter
    def partitions(self, value):
        self._partitions = value

    @property
    def dict(self):
        d = dict((item.lstrip('_'), getattr(self, item)) for item in self.__slots__)
        if d['partitions'] is not None:
            d['partitions'] = {i: subpar.dict for (i, subpar) in enumerate(d['partitions'])}
        return d

    def construct_from_dict(self, dict):
        super(Parameters, self).construct_from_dict(dict)
        try:
            tmp = {int(k): v for (k, v) in self._partitions.items()}
        except AttributeError:
            tmp = {}
        self._partitions = []
        for k in sorted(tmp):
            pp = PartitionParameters()
            pp.construct_from_dict(tmp[k])
            self.partitions.append(pp)
