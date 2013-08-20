#!/usr/bin/env python

from collection import Collection
from optimiser import Optimiser
import tempfile
import string
import random
import csv

import argparse

parser = argparse.ArgumentParser(description='Clustering optimiser')
parser.add_argument('-n', '--nclusters', type=int)
parser.add_argument('-f', '--format', default='phylip')
parser.add_argument('-d', '--datatype', default='protein')
parser.add_argument('-i', '--input_dir', default='./')
parser.add_argument('-c', '--compression', default=None)
parser.add_argument('-t', '--tmpdir', default='/tmp/')
# Collect all args for optimse and parse them later?
parser.add_argument('-r', '--nreassign', default=10, type=int)
parser.add_argument('-s', '--sample_size', default=10, type=int)
parser.add_argument('-o', '--output', default=None)

args = parser.parse_args()


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))

new_tmpdir = tempfile.mkdtemp(prefix='tmpwrap_mgp_', dir=args.tmpdir)

c = Collection(input_dir=args.input_dir,
               compression=args.compression, file_format=args.format, datatype=args.datatype,
               tmpdir=new_tmpdir)

o = Optimiser(args.nclusters, c)
o.optimise(max_iter=500, nreassign=args.nreassign, sample_size=args.sample_size)

output_name = args.output or 'output_' + id_generator(6)
output_fh = open(output_name, 'w+')
headings = ['Iteration', 'CPU Time', 'Likelihood', 'Partition']
output = [[i] + x for i, x in enumerate(o.Scorer.history)]

writer = csv.writer(output_fh, delimiter='\t', quoting=csv.QUOTE_NONE)
writer.writerow(headings)
writer.writerows(output)
