#!/usr/bin/env python
import sys
sys.path.append('/home/malcolm/Documents/EBi/Git/')
from treeCl.collection import Collection, Scorer
from treeCl.clustering import Clustering
from lib.remote.errors import optioncheck

import random
import tempfile
import string
import csv

import argparse

parser = argparse.ArgumentParser(description='Clustering optimiser')
parser.add_argument('-n', '--nclusters', type=int)
parser.add_argument('-f', '--format', default='phylip')
parser.add_argument('-d', '--datatype', default='protein')
parser.add_argument('-i', '--input_dir', default='./')
parser.add_argument('-c', '--compression', default=None)
parser.add_argument('-t', '--tmpdir', default='/tmp/')
parser.add_argument('-m', '--method', default='s')
parser.add_argument('-o', '--output', default=None)

args = parser.parse_args()

optioncheck(args.method, ['s', 'spectral', 'h', 'hierarchical', 'k', 'kmedoids', 'MDS', 'mds'])


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

new_tmpdir = tempfile.mkdtemp(prefix='tmpwrap_mgp_', dir=args.tmpdir)

c = Collection(input_dir=args.input_dir,
               compression=args.compression, file_format=args.format, datatype=args.datatype,
               tmpdir=new_tmpdir)

c.calc_NJ_trees()

dm = c.distance_matrix('euc')
cl = Clustering(dm)

if args.method in ['spectral', 's']:
    print 'Starting spectral decomposition...'
    decomp = cl.spectral_decomp('estimate', 'estimate')
    print 'Starting spectral clustering...'
    p = cl.spectral_cluster(4, decomp)
elif args.method in ['MDS', 'mds']:
    print 'Starting MDS decomposition...'
    decomp = cl.MDS_decomp()
    print 'Starting MDS clustering...'
    p = cl.MDS_cluster(4, decomp)
elif args.method in ['hierarchical', 'h']:
    print 'Starting hierarchical clustering...'
    p = cl.hierarchical(4, 'ward')
elif args.method in ['kmedoids', 'k']:
    print 'why would you do this'
    exit(0)

sc = Scorer(c.records, c.analysis)

score = sc.score(p)

output_name = args.output or 'output_' + id_generator(6)
output_fh = open(output_name, 'w+')
headings = ['Likelihood', 'Partition']
output = [score, p.partition_vector]

writer = csv.writer(output_fh, delimiter='\t', quoting=csv.QUOTE_NONE)
writer.writerow(headings)
writer.writerow(output)
