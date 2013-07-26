#!/usr/bin/env python

from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import os


PATTERNS = {'uniq': re.compile(r'[A-Z0-9]{6}$')}


class Result(object):
    def __init__(self, likelihood, clusters, history, uniq_id=None, info=None):
        self.likelihood = likelihood
        self.clusters = clusters
        self.cputime = history[-1][0]
        self.history = history
        self.id = uniq_id
        self.info = info

    def __len__(self):
        return(len(self.history))

    def timeseries(self):
        return([(round(rec[0], 1), int(rec[1])) for rec in self.history])

    def by_iteration(self):
        return([(i, int(rec[1])) for i, rec in enumerate(self.history)])

    def print_table(self, file=sys.stdout):
        for i, rec in enumerate(self.history):
            print('i\t', 'cpu time\t', 'likelihood')
            print(i, round(rec[0], 1), int(rec[1]), sep='\t', end='\n')

    # @property
    # def plot(self):
    #     return(self._plot)

    # @plot.setter
    def plot(self, name, ext='png', time=True):
        if time:
            data = np.asarray(self.timeseries())
        else:
            data = np.asarray(self.by_iteration())
        plt.plot(data[:, 0], data[:, 1])  # WHERE IS PLOT OBJECT STORED
        filename = name + '.' + ext
        plt.xlabel = 'System Time'
        plt.ylabel = 'Log Likelihood'
        plt.savefig(filename)
        # self._plot = plt.plot(data[:, 0], data[:, 1], , ylab='Log Likelihood of Partition')


def fileparser(f):
    uniq_id = PATTERNS['uniq'].search(f.name).group(0)
    clusters = []
    history = []
    _ = f.readline()
    best_score = -np.Inf

    while True:
        line = f.readline()
        try:
            niter, time, score, clusters = [eval(x) for x in line.rstrip().split('\t')]
            if score > best_score:
                best_score = score
                final_partition = clusters
        except:
            break
        history.append((time, score))

    result = Result(best_score, final_partition, history, uniq_id)
    return(result)


def foldersearch(folderpath, info=None, prefix='output_'):
    filenames = os.listdir(folderpath)
    results = []
    for fn in filenames:
        if not re.match(prefix, fn):
            continue
        f = open(os.path.join(folderpath, fn))
        result = fileparser(f)
        if info:
            result.info = info
        results.append(result)

    return(results)


def plot_folder(folder):
    results = foldersearch(folder)

    plt.hold(False)

    for r in results:
        plotname = os.path.join(folder, 'plot_' + str(r.id))
        r.plot(plotname, 'png')

    lls = [r.likelihood for r in results]
    times = [r.cputime for r in results]

    print('For Folder: {0}'.format(folder))
    print('Mean Likelihood:' + str(mean(lls)))
    print('Mean time:' + str(mean(times)))

    plt.scatter(times, lls)
    plt.savefig(os.path.join(folder, 'scatterplot'))

if __name__ == '__main__':
    from numpy import mean
    print('Running as script.')

    try:
        folders = sys.argv[1:]
    except:
        folder = './'

    print(folders)

    for folder in folders:
        plot_folder(folder)
