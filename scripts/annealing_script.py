#!/usr/bin/env python

from collection import Collection, Scorer
from clustering import Partition
from random import randint
from anneal import *
import pickle

c = Collection(input_dir='/home/malcolm/Documents/EBi/Data/easy_case', compression='gz', file_format='phylip', datatype='protein')

scorer = Scorer(c.records, 'nj')

k = 4

partition = [randint(1, k) for rec in scorer.records]


def likelihood(partition, scorer):
    score = scorer.score(Partition(partition))
    return(score)

print type(partition)

opts = {'func': likelihood,
        'x0': partition,
        'args': [scorer],
        'schedule': 'cluster',
        'full_output': 1,
        'T0': 100000,
        'Tf': 1,
        'maxeval': None,
        'maxaccept': None,
        'maxiter': 400,
        'boltzmann': 1.0,
        'learn_rate': 0.5,
        'feps': 1e-6,
        'quench': 1.0,
        'm': 1.0,
        'n': 1.0,
        'lower': -100,
        'upper': 100,
        'dwell': 50,
        'disp': True,
        }

res = anneal(**opts)

pickle.dumps(res, 'SA_results')
history = scorer.history()
pickle.dumps(history, 'SA_history')