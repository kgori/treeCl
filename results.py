#!/usr/bin/env python

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import os


class Result(object):
    def __init__(self, likelihood, clusters, history):
        self.likelihood = likelihood
        self.clusters = clusters
        self.history = history

    def __len__(self):
        return(len(self.history))

    def timeseries(self):
        return(self.history)

    def by_iteration(self):
        return([(i, rec[1]) for i, rec in enumerate(self.history)])

    def table(self, file=sys.stdout):
        for i, rec in enumerate(self.history):
            print(i, rec[0], rec[1], sep='\t', end='\n')

    @property
    def plot(self):
        return(self._plot)

    @plot.setter
    def plot(self, time=True):
        if time:
            data = np.asarray(self.timeseries())
        else:
            data = np.asarray(self.by_iteration())
        self._plot = plt.plot(data[:, 0], data[:, 1], xlab='System Time', ylab='Log Likelihood of Partition')

    def print_plot(self, name, type='png'):
        filename = name + '.' + extension
        pass
        # Need to know how to control muliple plots, where plot objects are stored etc - need to copy plot???


def fileparser(f):
    likelihood = eval(re.match(r'Likelihood score:(\-?\d+\.\d+)', f.readline()).group(1))

    clusters = []
    while True:
        line = f.readline()
        try:
            cluster = re.search(r'(\(.*\))', line).group(1)
        except:
            break
        clusters.append(eval(cluster))

    history = []
    while True:
        line = f.readline()
        try:
            niter, time, score = [eval(x) for x in line.rstrip().split('\t')]
        except:
            break
        history.append((time, score))

    result = Result(likelihood, clusters, history)
    return(result)


def foldersearch(folder):
    files = os.listdir(folder)
    results = []
    for f in files:
        res = open(os.path.join(folder, f))
        results.append(fileparser(res))

    return(results)
