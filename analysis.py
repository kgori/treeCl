#!/usr/bin/env python

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import os


PATTERNS = {'ll': re.compile(r'Likelihood score:(\-?\d+\.\d+)'),
            'cluster': re.compile(r'(\(.*\))'),
            'uniq': re.compile(r'[A-Z0-9]{6}$')}


class Result(object):
    def __init__(self, likelihood, clusters, history, uniq_id=None, info=None):
        self.likelihood = likelihood
        self.clusters = clusters
        self.history = history
        self.id = uniq_id
        self.info = info

    def __len__(self):
        return(len(self.history))

    def timeseries(self):
        return(self.history)

    def by_iteration(self):
        return([(i, rec[1]) for i, rec in enumerate(self.history)])

    def print_table(self, file=sys.stdout):
        for i, rec in enumerate(self.history):
            print(i, rec[0], rec[1], sep='\t', end='\n')

    # @property
    # def plot(self):
    #     return(self._plot)

    # @plot.setter
    def plot(self, name, ext='png', time=True):
        if time:
            data = np.asarray(self.timeseries())
        else:
            data = np.asarray(self.by_iteration())
        plt.plot(data[:, 0], data[:, 1]) # WHERE IS PLOT OBJECT STORED
        filename = name + '.' + ext
        plt.xlabel='System Time'
        plt.ylabel='Log Likelihood'
        plt.savefig(filename)
        # self._plot = plt.plot(data[:, 0], data[:, 1], , ylab='Log Likelihood of Partition')


def fileparser(f):
    print(f.name)
    uniq_id = PATTERNS['uniq'].search(f.name).group(0)
    likelihood = eval(PATTERNS['ll'].match(f.readline()).group(1))
    clusters = []
    history = []

    while True:
        line = f.readline()
        try:
            cluster = PATTERNS['cluster'].search(line).group(1)
        except:
            break
        clusters.append(eval(cluster))

    while True:
        line = f.readline()
        try:
            niter, time, score = [eval(x) for x in line.rstrip().split('\t')]
        except:
            break
        history.append((time, score))

    result = Result(likelihood, clusters, history, uniq_id)
    return(result)


def foldersearch(folderpath, info=None):
    filenames = os.listdir(folderpath)
    results = []
    print(filenames)
    for fn in filenames:
        f = open(os.path.join(folderpath, fn))
        result = fileparser(f)
        if info:
            result.info = info
        results.append(result)

    return(results)
