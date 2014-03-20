#!/usr/bin/env python
from __future__ import print_function
import tree_collection

mat = open('/Users/kgori/git/treesoflife/data/DistVar.txt').read()
map_= open('/Users/kgori/git/treesoflife/data/GenomeMap.txt').read()
lab = open('/Users/kgori/git/treesoflife/data/Labels.txt').read()
tre = open('/Users/kgori/git/treesoflife/data/Tree.nwk').read()

print('Running...')
result = tree_collection.compute(mat, map_, lab, tre, 3, False)
print(type(result))
print(result[0])
print(result[1])