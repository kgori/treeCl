# from __future__ import print_function

from alignment import Alignment
from clustering import Spectral, Hierarchical, MultidimensionalScaling, Automatic, Evaluation
from collection import Collection, Scorer
from concatenation import Concatenation
from distance_matrix import DistanceMatrix
from optimiser import Optimiser, EM
from partition import Partition
from plotter import Plotter
from simulator import Simulator
from tree import Tree

import logging.config
import yaml
from pkg_resources import resource_string
conf = resource_string(__name__, 'logging/logging.yaml')

D = yaml.load(conf)
D.setdefault('version', 1)
logging.config.dictConfig(D)
del D
