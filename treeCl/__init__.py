# from __future__ import print_function

from alignment import Alignment
from clustering import Clustering, Spectral
from collection import Collection, Scorer
from concatenation import Concatenation
from distance_matrix import DistanceMatrix
from optimiser import Optimiser, EM
from partition import Partition
from plotter import Plotter
from simulator import Simulator
from tree import Tree

import logging
import yaml
from pkg_resources import resource_string
logging_config = resource_string(__name__, 'logging_config')




D = yaml.load(logging_config)
D.setdefault('version', 1)
logging.config.dictConfig(D)
logger = logging.getLogger(__name__)
logger.info('xyz')