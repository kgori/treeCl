#!/usr/bin/env python

import analysis
import matplotlib.pyplot as plt
from numpy import mean
import sys
import os

try:
    folders = sys.argv[1:]
except:
    folder = './'


def plot(folder):
    results = analysis.foldersearch(folder)

    plt.hold(False)

    for r in results:
        plotname = os.path.join(folder, 'plot_' + str(r.id))
        r.plot(plotname, 'png')

    lls = [r.likelihood for r in results]
    times = [r.cputime for r in results]

    print 'For folder {}:'.format(folder)
    print 'Mean Likelihood:' + str(mean(lls))
    print 'Mean time:' + str(mean(times))

    plt.scatter(times, lls)
    plt.savefig(os.path.join(folder, 'scatterplot'))

for folder in folders:
    plot(folder)
