#!/usr/bin/env python


from ..externals.gtp import GTP
from ...remote.datastructs.tree import Tree

class TrClTree(Tree):

    def geodist(self, other):
        gtp = GTP()
        return gtp.pairwise(self, other)

    @classmethod
    def cast(cls, tree):
        """ Copies Tree object as TrClTree object """
        cast = TrClTree().clone_from(tree)
        return cast

    @classmethod
    def new_yule(
        cls,
        nspecies,
        names=None,
        cf=False,
        ):
        t = super(TrClTree, cls).new_yule(nspecies, names, cf)
        return cls.cast(t)

    @classmethod
    def new_coal(
        cls,
        nspecies,
        names=None,
        cf=False,
        ):
        t = super(TrClTree, cls).new_coal(nspecies, names, cf)
        return cls.cast(t)

    @classmethod
    def new_rtree(
        cls,
        nspecies,
        names=None,
        cf=False,
        ):
        t = super(TrClTree, cls).new_rtree(nspecies, names, cf)
        return cls.cast(t)
