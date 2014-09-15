from __future__ import print_function
from collection import Collection, Scorer, Concatenation
from clustering import Clustering, Partition
from distance_matrix import DistanceMatrix
from optimiser import Optimiser
from plotter import Plotter
from simulator import Simulator
from tree import Tree
from utils import fileIO

def set_java_classpath():
    import os
    current_classpath = os.getenv('CLASSPATH')
    jarfile_location = os.path.join(os.path.dirname(__file__), 'software_interfaces', 'gtp.jar')
    if current_classpath is not None:
        if 'gtp.jar' in current_classpath:
            return
        else:
            os.environ['CLASSPATH'] = ':'.join([current_classpath, jarfile_location])
    else:
        os.environ['CLASSPATH'] = jarfile_location

def setup_java_classes():
    set_java_classpath()
    from jnius import autoclass
    PhyloTree = autoclass('distanceAlg1.PhyloTree')
    GTP = autoclass('polyAlg.PolyMain')
    return PhyloTree, GTP

PhyloTree, GTP = setup_java_classes()
