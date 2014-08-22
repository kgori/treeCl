#!/usr/bin/env python

# treeCl
from tree import Tree
from ..constants import TMPDIR
from ..software_interfaces.gtp import GTP


class TrClTree(Tree):
    def geodist(self, other, tmpdir=None, normalise=False):
        tmpdir = tmpdir or TMPDIR
        gtp = GTP(tmpdir=tmpdir)
        if self ^ other:
            t1, t2 = self.__class__.pruned_pair(self, other)
        else:
            t1, t2 = self, other
        if normalise:
            t1, t2 = self.__class__.normalised_pair(t1, t2)
        return gtp.pairwise(t1, t2)

    @classmethod
    def cast(cls, tree):
        """ Copies Tree object as TrClTree object """
        cast = TrClTree().clone_from(tree)
        return cast

    @classmethod
    def new_iterative_rtree(cls, nspecies):
        t = super(TrClTree, cls).new_iterative_rtree(nspecies)
        return cls.cast(t)

    @classmethod
    def new_yule(cls, nspecies, **kwargs):
        t = super(TrClTree, cls).new_yule(nspecies, **kwargs)
        return cls.cast(t)

    @classmethod
    def new_coal(cls, nspecies, **kwargs):
        t = super(TrClTree, cls).new_coal(nspecies, **kwargs)
        return cls.cast(t)

    @classmethod
    def new_rtree(cls, nspecies, **kwargs):
        t = super(TrClTree, cls).new_rtree(nspecies, **kwargs)
        return cls.cast(t)
