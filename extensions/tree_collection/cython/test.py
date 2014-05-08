#!/usr/bin/env python
from __future__ import print_function
import tree_collection

mat = open('/Users/kgori/git_repos/gitolite/treesoflife/data/DistVar.txt').read()
map_= open('/Users/kgori/git_repos/gitolite/treesoflife/data/GenomeMap.txt').read()
lab = open('/Users/kgori/git_repos/gitolite/treesoflife/data/Labels.txt').read()
tre = open('/Users/kgori/git_repos/gitolite/treesoflife/data/Tree.nwk').read()

print('Running test 1 - fit initial tree...')
result = tree_collection.compute(mat, map_, lab, tre,
                                 iter=8,
                                 keep_topology=True,
                                 quiet=False)
print(type(result))
print(result[0])
print(result[1])

print('Running test 2 - optimise tree...')
result = tree_collection.compute(mat, map_, lab, tre,
                                 iter=8,
                                 keep_topology=False,
                                 quiet=False)
print(type(result))
print(result[0])
print(result[1])
