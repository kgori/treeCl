#!/usr/bin/env python
from __future__ import print_function

# standard library
import itertools
import random
import re

# third party
import dendropy
import numpy as np
from tree_distance import PhyloTree

# treeCl
from errors import optioncheck
from utils import fileIO
from utils.decorators import lazyprop



def cast(dendropy_tree):
    """ Cast dendropy.Tree instance as Tree instance """
    return Tree(dendropy_tree.as_newick_string() + ';')


def _infinite_labels_generator(labels, start=2, shuffle=True):
    l = len(labels)
    loop1 = random.sample(labels, l) if shuffle else labels
    return itertools.chain.from_iterable([loop1, ('{}{}'.format(x, y) for x, y in
                                                  itertools.izip(itertools.cycle(labels),
                                                                 itertools.chain.from_iterable(
                                                                     itertools.repeat(i, len(loop1)) for i
                                                                     in
                                                                     itertools.count(start, 1))))])


def edge_length_check(length, edge):
    """ Raises error if length is not in interval [0, edge.length] """
    try:
        assert 0 <= length <= edge.length
    except AssertionError:
        if length < 0:
            raise TreeError('Negative edge-lengths are disallowed')
        raise TreeError(
            'This edge isn\'t long enough to prune at length {0}\n'
            '(Edge length = {1})'.format(length, edge.length))


def rootcheck(edge, msg='This is the root edge'):
    """ Raises error if edge is the root edge (has no tail node) """
    if not edge.tail_node:
        raise TreeError(msg)


def logn_correlated_rate(parent_rate, branch_length, autocorrel_param, size=1):
    """
    The log of the descendent rate, ln(Rd), is ~ N(mu, bl*ac), where
    the variance = bl*ac = branch_length * autocorrel_param, and mu is set
    so that E[Rd] = Rp:
    E[X] where ln(X) ~ N(mu, sigma^2) = exp(mu+(1/2)*sigma_sq)
    so Rp = exp(mu+(1/2)*bl*ac),
    ln(Rp) = mu + (1/2)*bl*ac,
    ln(Rp) - (1/2)*bl*ac = mu,
    so ln(Rd) ~ N(ln(Rp) - (1/2)*bl*ac, bl*ac)
    (NB: Var[Rd] = Rp^2 * (exp(bl*ac)-1),
         Std[Rd] = Rp * sqrt(exp(bl*ac)-1)

    See: H Kishino, J L Thorne, and W J Bruno (2001)
    """
    if autocorrel_param <= 0:
        raise Exception('Autocorrelation parameter must be greater than 0')

    variance = branch_length * autocorrel_param
    stdev = np.sqrt(variance)
    ln_descendant_rate = np.random.normal(np.log(parent_rate) - 0.5 * variance,
                                          scale=stdev, size=size)
    descendant_rate = np.exp(ln_descendant_rate)
    return float(descendant_rate) if size == 1 else descendant_rate


class TreeError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return self.msg


class SPR(object):
    """ Subtree prune and regraft functionality """

    def __init__(self, tree):
        self.tree = tree

    def _check_single_outgroup(self):
        """
        If only one (or none) of the seed node children is not a leaf node
        it is not possible to prune that edge and make a topology-changing
        regraft.
        """
        root_child_nodes = self.tree.seed_node.child_nodes()
        not_leaves = np.logical_not([n.is_leaf() for n in root_child_nodes])
        if not_leaves[not_leaves].size <= 1:
            return [root_child_nodes[np.where(not_leaves)[0]].edge]
        return []

    def prune(self, edge, length=None):
        """ Prunes a subtree from the main Tree, retaining an edge length
        specified by length (defaults to entire length). The length is sanity-
        checked by edge_length_check, to ensure it is within the bounds
        [0, edge.length].

        Returns the basal node of the pruned subtree. """

        length = length or edge.length
        edge_length_check(length, edge)

        n = edge.head_node
        self.tree.prune_subtree(n, delete_outdegree_one=False)
        n.edge_length = length
        self.tree._dirty = True
        return n

    def regraft(self, edge, node, length=None):
        """ Grafts a node onto an edge of the Tree, at a point specified by
        length (defaults to middle of edge). """

        rootcheck(edge, 'SPR regraft is not allowed on the root edge')
        length = length or edge.length / 2.  # Length measured from head to tail
        edge_length_check(length, edge)

        t = edge.tail_node
        h = edge.head_node
        new = t.new_child(edge_length=edge.length - length)
        t.remove_child(h)
        new.add_child(h, edge_length=length)
        new.add_child(node)
        self.tree._dirty = True
        self.tree.update_splits(delete_outdegree_one=True)

    def spr(self, prune_edge, length1, regraft_edge, length2):
        assert (regraft_edge.head_node
                not in prune_edge.head_node.preorder_iter())
        node = self.prune(prune_edge, length1)
        self.regraft(regraft_edge, node, length2)
        self.tree._dirty = True

    def rspr(self, disallow_sibling_sprs=False,
             keep_entire_edge=False, rescale=False):
        """ Random SPR, with prune and regraft edges chosen randomly, and
        lengths drawn uniformly from the available edge lengths.

        N1: disallow_sibling_sprs prevents sprs that don't alter the topology
        of the tree """

        starting_length = self.tree.length()

        excl = [self.tree.seed_node.edge]  # exclude r
        if disallow_sibling_sprs:
            excl.extend(self._check_single_outgroup())
        prune_edge, l1 = self.tree.map_event_onto_tree(excl)
        if keep_entire_edge:
            l1 = prune_edge.length

        prune_edge_child_nodes = prune_edge.head_node.preorder_iter()
        excl.extend([node.edge for node in prune_edge_child_nodes])

        if disallow_sibling_sprs:
            sibs = [node.edge for node in prune_edge.head_node.sister_nodes()]
            par = prune_edge.tail_node.edge
            sibs.append(par)
            for edge in sibs:
                if edge not in excl:
                    excl.append(edge)
            if set(self.tree.preorder_edge_iter()) - set(excl) == set([]):
                print(repr(self.tree))
                print(self.tree.as_ascii_plot())
                # print(edges[prune_edge])
                raise Exception('No non-sibling sprs available')

        regraft_edge, l2 = self.tree.map_event_onto_tree(excl)

        self.spr(prune_edge, l1, regraft_edge, l2)
        if rescale:
            self.tree.scale(starting_length / self.tree.length())
            self.tree._dirty = True


class LGT(object):
    def __init__(self, tree):
        self.SPR = SPR(tree)
        self.tree = self.SPR.tree
        try:
            self.tree.calc_node_ages()
        except:
            raise Exception('Tree is not ultrametric')

    def get_time(self, *args):
        e, l = self.tree.map_event_onto_tree(*args)
        time = l + e.head_node.age
        return time

    def matching_edges(self, time):
        def edge_matches_time(edge):
            if edge.tail_node is None:
                return False
            return edge.head_node.age < time < edge.tail_node.age

        matching_edges = self.tree.preorder_edge_iter(edge_matches_time)
        return list(matching_edges)

    def rlgt(self, time=None, disallow_sibling_lgts=False):
        self.tree.calc_node_ages()
        excl = [self.tree.seed_node.edge]

        if time is None:
            if disallow_sibling_lgts:
                self.add_single_node()
                children = self.tree.seed_node.child_nodes()
                excl.extend([n.edge for n in children])
                time = self.get_time(excl)
                print('time = {0}'.format(time))
                self.tree.update_splits()
                self.tree._dirty = True

            else:
                time = self.get_time(excl)
                print('time = {0}'.format(time))

        matching_edges = self.matching_edges(time)
        donor = random.sample(matching_edges, 1)[0]
        matching_edges.remove(donor)
        if disallow_sibling_lgts:
            sibs = donor.head_node.sister_nodes()
            for sib in sibs:
                if sib.edge in matching_edges:
                    matching_edges.remove(sib.edge)
        receiver = random.sample(matching_edges, 1)[0]

        l1 = time - receiver.head_node.age
        l2 = time - donor.head_node.age

        self.SPR.spr(receiver, l1, donor, l2)
        self.tree.calc_node_ages()

    def add_single_node(self):
        cn = self.tree.seed_node.child_nodes()
        el = lambda n: n.edge_length
        sh = min(cn, key=el)
        lo = max(cn, key=el)
        new = self.tree.seed_node.new_child(edge_length=sh.edge_length)
        self.tree.prune_subtree(lo, delete_outdegree_one=False)
        lo.edge_length -= sh.edge_length
        new.add_child(lo)
        self.tree.update_splits(delete_outdegree_one=False)
        self.tree.calc_node_ages()
        self.tree._dirty = True


class NNI(object):
    def __init__(self, tree):
        self.tree = tree
        if tree.rooted:
            self.reroot = True
            self.rooting_info = self.tree.reversible_deroot()
        else:
            self.reroot = False
            self.rooting_info = None

    def get_children(self, inner_edge):
        """ Given an edge in the tree, returns the child nodes of the head and
        the tail nodes of the edge, for instance:

            A      C    | A, B, C and D are the children of the edge --->,
             \    /     | C and D are the head node children, and A and B
             t--->h     | are the tail node children.
             /    \
            B      D    | Output: {'head': [<C>, <D>], 'tail': [<A>, <B>]}

        N1: Edges are directional in dendropy trees. The head node of an
        edge is automatically a child of the tail node, but we don't want this.
        """

        h = inner_edge.head_node
        t = inner_edge.tail_node
        if not self.tree.seed_node == t:
            original_seed = self.tree.seed_node
            self.tree.reseed_at(t)
        else:
            original_seed = None
        head_children = h.child_nodes()
        tail_children = list(set(t.child_nodes()) - {h})  # See N1
        if original_seed:
            self.tree.reseed_at(original_seed)

        return {'head': head_children, 'tail': tail_children}

    def nni(
            self,
            edge,
            head_subtree,
            tail_subtree,
    ):
        """ *Inplace* Nearest-neighbour interchange (NNI) operation.

        An edge in the tree has two or more subtrees at each end (ends are
        designated 'head' and 'tail'). The NNI operation exchanges one of the
        head subtrees for one of the tail subtrees, as follows:

            A      C                        C      A    | Subtree A is exchanged
             \    /        +NNI(A,C)         \    /     | with subtree C.
              --->        ==========>         --->      |
             /    \                          /    \     |
            B      D                        B      D


        """

        # This implementation works on unrooted Trees. If the input Tree is
        # rooted, the ReversibleDeroot decorator will temporarily unroot the
        # tree while the NNI is carried out

        original_seed = self.tree.seed_node
        head = edge.head_node
        tail = edge.tail_node
        self.tree.reseed_at(tail)
        try:
            assert head_subtree.parent_node == head
            assert tail_subtree.parent_node == tail
        except:
            print(head, tail, head_subtree, tail_subtree)
            raise
        head.remove_child(head_subtree)
        tail.remove_child(tail_subtree)
        head.add_child(tail_subtree)
        tail.add_child(head_subtree)
        self.tree.reseed_at(original_seed)
        self.tree.update_splits()
        self.tree._dirty = True

    def reroot_tree(self):
        if self.reroot and self.rooting_info is not None:
            self.tree.reroot_at_edge(*self.rooting_info)
            self.tree.update_splits()
            self.tree._dirty = True
        return self.tree

    def rnni(self):
        e = random.choice(self.tree.get_inner_edges())
        children = self.get_children(e)
        h = random.choice(children['head'])
        t = random.choice(children['tail'])
        self.nni(e, h, t)


class Tree(dendropy.Tree):
    """Augmented version of dendropy Tree class"""

    def __init__(
            self,
            newick=None,
            name=None,
            **kwargs
    ):

        super(Tree, self).__init__(**kwargs)
        if newick:
            self.read_from_string(newick, 'newick', **kwargs)
            if self.rooted:
                self.is_rooted = True
                self.update_splits()
        self.name = name
        self._phylotree = None
        self._dirty = False

    def __repr__(self):
        return '{0}{1}'.format(self.__class__.__name__,
                               (self.newick if self.newick else '(None)'))

    def __str__(self):
        """ Represents the object's information inside a newick comment, so is
        still interpretable by a (good) newick parser """

        s = 'Tree Object: {}\n'.format(self.name)
        s += self.newick

        return s

    def __len__(self):
        """ Number of leaves on the Tree. For total branch length use
        self.length()"""
        return len(self.leaf_nodes())

    def __and__(self, other):
        """ Overloads & operator:

        'self & other' is equivalent to 'self.intersection(other)''
        """
        return self.intersection(other)

    def __xor__(self, other):
        return self.labels ^ other.labels

    @property
    def labels(self):
        """ Returns the taxon set of the tree (same as the label- or
        leaf-set) """
        return set([n.taxon.label for n in self.leaf_nodes()])

    def sample_labels(self, n):
        """ Returns a set of n labels sampled from the labels of the tree
        :param n: Number of labels to sample
        :return: set of randomly sampled labels
        """
        if n >= len(self):
            return self.labels
        sample = random.sample(self.labels, n)
        return set(sample)

    @property
    def newick(self):
        """
        For more control the dendropy method self.as_string('newick', **kwargs)
        can be used.
        KWargs include:
            internal_labels [True/False]  - turn on/off bootstrap labels
            suppress_rooting [True/False] - turn on/off [&U] or [&R] rooting
                                            state labels
            edge_label_compose_func      - function to convert edge lengths:
                                            takes edge as arg, returns string
        """
        n = self.as_newick_string()
        if n:
            return n if n.endswith(';') else n + ';'
        return n

    @property
    def phylotree(self):
        """
        Gets the c++ PhyloTree object corresponding to this tree.
        Should be canonically the same - we set a _dirty flag if the Python version of the
        tree has changed since construction. If the flag is set then we reconstruct
        the c++ PhyloTree
        :return: PhyloTree instance
        """
        if not self._phylotree or self._dirty:
            self._phylotree = PhyloTree(self.newick, self.rooted)
            self._dirty = False
        return self._phylotree

    @newick.setter
    def newick(self, newick_string):
        if self.newick:
            print('Newick string already loaded: {0}'.format(self.newick))
            return
        self.read_from_string(newick_string, 'newick')

    @property
    def rooted(self):
        """ Predicate testing for rootedness by checking for a bifurcation
        at the root. """
        return len(self.seed_node.child_nodes()) == 2 if self.newick else None

    @classmethod
    def bifurcate_base(cls, newick):
        """ Rewrites a newick string so that the base is a bifurcation
        (rooted tree) """
        t = cls(newick)
        t.resolve_polytomies()
        return t.newick

    @classmethod
    def trifurcate_base(cls, newick):
        """ Rewrites a newick string so that the base is a trifurcation
        (usually means an unrooted tree) """
        t = cls(newick)
        t.deroot()
        return t.newick

    def copy(self):
        """ Returns an independent copy of self """
        return self.__class__(self.newick)

    def get_inner_edges(self):
        """ Returns a list of the internal edges of the tree. """
        inner_edges = [e for e in self.preorder_edge_iter() if e.is_internal()
                       and e.head_node and e.tail_node]
        return inner_edges

    def get_nonroot_edges(self):
        return [e for e in self.preorder_edge_iter()
                if e.head_node and e.tail_node]

    def intersection(self, other):
        """ Returns the intersection of the taxon sets of two Trees """
        taxa1 = self.labels
        taxa2 = other.labels
        return taxa1 & taxa2

    def map_event_onto_tree(self, excluded_edges=None):
        edge_list = list(self.preorder_edge_iter())
        if excluded_edges is not None:
            if not isinstance(excluded_edges, list):
                excluded_edges = [excluded_edges]
            for excl in excluded_edges:
                try:
                    edge_list.remove(excl)
                except ValueError:
                    print('Excluded_edges list includes some things')
                    print('that aren\'t in the tree')
                    print('like this one:', excl)
        lengths = np.array([edge.length for edge in edge_list])
        cumulative_lengths = lengths.cumsum()
        rnum = np.random.random() * cumulative_lengths[-1]
        index = cumulative_lengths.searchsorted(rnum)
        chosen_edge = edge_list[index]
        from_head_length = cumulative_lengths[index] - rnum
        return chosen_edge, from_head_length

    def multifurcate(self, threshold=1e-06, update_splits=True):
        for edge in self.postorder_edge_iter():
            if edge.is_internal():
                if edge.length <= threshold:
                    edge.collapse()
                    self._dirty = True
        if update_splits:
            self.update_splits()

    def ntaxa(self):
        return len(self)

    def pairdist(self, taxon_label1, taxon_label2):
        if self.patristic is None:
            print('Error calculating patristic distances - maybe this '
                  'tree has no branch lengths?')
            return

        leaf1 = self.find_node_with_taxon_label(taxon_label1)
        leaf2 = self.find_node_with_taxon_label(taxon_label2)

        if leaf1:
            taxon1 = leaf1.taxon
        else:
            print('Couldn\'t find {0} on the tree'.format(taxon_label1))
            return

        if leaf2:
            taxon2 = leaf2.taxon
        else:
            print('Couldn\'t find {0} on the tree'.format(taxon_label2))
            return

        return self.patristic(taxon1, taxon2)

    @lazyprop
    def patristic(self):
        try:
            pdm = dendropy.treecalc.PatristicDistanceMatrix(self)
        except TypeError:
            pdm = None
        return pdm

    def prune_to_subset(self, subset, inplace=False):
        """ Prunes the Tree to just the taxon set given in `subset` """
        if not subset.issubset(self.labels):
            print('"subset" is not a subset')
            return
        if not inplace:
            t = self.copy()
        else:
            t = self
        t.retain_taxa_with_labels(subset)
        t.update_splits()
        t._dirty = True
        return t

    def randomise_branch_lengths(
            self,
            i=(1, 1),
            l=(1, 1),
            distribution_func=random.gammavariate,
            inplace=False,
    ):
        """ Replaces branch lengths with values drawn from the specified
        distribution_func. Parameters of the distribution are given in the
        tuples i and l, for interior and leaf nodes respectively. """

        if not inplace:
            t = self.copy()
        else:
            t = self

        for n in t.preorder_node_iter():
            if n.is_internal():
                n.edge.length = max(0, distribution_func(*i))
            else:
                n.edge.length = max(0, distribution_func(*l))
        t._dirty = True
        return t

    def randomise_labels(
            self,
            inplace=False,
    ):
        """ Shuffles the leaf labels, but doesn't alter the tree structure """

        if not inplace:
            t = self.copy()
        else:
            t = self

        names = t.labels
        random.shuffle(list(names))
        for l in t.leaf_iter():
            l.taxon_label = names.pop()
        t._dirty = True
        return t

    def reversible_deroot(self):
        """ Stores info required to restore rootedness to derooted Tree. Returns
        the edge that was originally rooted, the length of e1, and the length
        of e2.

        Dendropy Derooting Process:
        In a rooted tree the root node is bifurcating. Derooting makes it
        trifurcating.

        Call the two edges leading out of the root node e1 and e2.
        Derooting with Tree.deroot() deletes one of e1 and e2 (let's say e2),
        and stretches the other to the sum of their lengths. Call this e3.

        Rooted tree:                   Derooted tree:
                 A                         A   B
                 |_ B                       \ /
                /                            |
               /e1                           |e3 (length = e1+e2; e2 is deleted)
        Root--o               ===>           |
               \e2                     Root--o _ C
                \ _ C                        |
                 |                           D
                 D

        Reverse this with Tree.reroot_at_edge(edge, length1, length2, ...)
        """
        root_edge = self.seed_node.edge
        lengths = dict([(edge, edge.length) for edge
                        in self.seed_node.incident_edges() if edge is not root_edge])
        self.deroot()
        reroot_edge = (set(self.seed_node.incident_edges())
                       & set(lengths.keys())).pop()
        self.update_splits()
        self._dirty = True
        return (reroot_edge, reroot_edge.length - lengths[reroot_edge],
                lengths[reroot_edge])

    def autocorrelated_relaxed_clock(self, root_rate, autocorrel,
                                     distribution='lognormal'):
        """
        Attaches rates to each node according to autocorrelated lognormal
        model from Kishino et al.(2001), or autocorrelated exponential
        """
        optioncheck(distribution, ['exponential', 'lognormal'])

        if autocorrel == 0:
            for node in self.preorder_node_iter():
                node.rate = root_rate
            return

        for node in self.preorder_node_iter():
            if node == self.seed_node:
                node.rate = root_rate
            else:
                parent_rate = node.parent_node.rate
                bl = node.edge_length
                if distribution == 'lognormal':
                    node.rate = logn_correlated_rate(parent_rate, bl,
                                                     autocorrel)
                else:
                    node.rate = np.random.exponential(parent_rate)

    def uncorrelated_relaxed_clock(self, root_rate, variance,
                                   distribution='lognormal'):
        optioncheck(distribution, ['exponential', 'lognormal'])

        for node in self.preorder_node_iter():
            if node == self.seed_node:
                node.rate = root_rate
            else:
                if distribution == 'lognormal':
                    mu = np.log(root_rate) - 0.5 * variance
                    node.rate = np.random.lognormal(mu, variance)
                else:
                    node.rate = np.random.exponential(root_rate)

    def rlgt(self, time=None, times=1,
             disallow_sibling_lgts=False):
        """ Uses class LGT to perform random lateral gene transfer on
        ultrametric tree """

        lgt = LGT(self.copy())
        for _ in range(times):
            lgt.rlgt(time, disallow_sibling_lgts)
        return lgt.tree

    def rnni(self, times=1):
        """ Applies a NNI operation on a randomly chosen edge. """

        nni = NNI(self.copy())
        for _ in range(times):
            nni.rnni()
        nni.reroot_tree()
        return nni.tree

    def rspr(self, times=1, **kwargs):
        """ Random SPR, with prune and regraft edges chosen randomly, and
        lengths drawn uniformly from the available edge lengths.

        N1: disallow_sibling_sprs prevents sprs that don't alter the topology
        of the tree """

        spr = SPR(self.copy())
        for _ in range(times):
            spr.rspr(**kwargs)
        return spr.tree

    def scale(self, factor, inplace=True):
        """ Multiplies all branch lengths by factor. """
        if not inplace:
            t = self.copy()
        else:
            t = self
        t.scale_edges(factor)
        t._dirty = True
        return t

    def strip(self, inplace=False):
        """ Sets all edge lengths to None """
        if not inplace:
            t = self.copy()
        else:
            t = self
        for e in t.preorder_edge_iter():
            e.length = None
        t._dirty = True
        return t

    @classmethod
    def read_from_file(cls, infile, taxon_set=None):
        with fileIO.freader(infile) as reader:
            s = reader.read()

        return cls.gen_from_text(s, taxon_set=taxon_set)

    # @classmethod
    # def gen_from_text(cls, s, taxon_set=None):
    # """ Generate new tree object from str method output """
    #
    #     name_search = re.search(r'(?<=Name:\t)(\w+)+', s)
    #     program_search = re.search(r'(?<=Program:\t)(\w+)+', s)
    #     score_search = re.search(r'(?<=Score:\t)([0-9.\-\+]+)', s)
    #     tree_search = re.search(r'(?<=Tree:\t]).+', s)
    #
    #     name = regex_search_extract(name_search)
    #     program = regex_search_extract(program_search)
    #     score = regex_search_extract(score_search)
    #     tree = regex_search_extract(tree_search)
    #
    #     if score:
    #         try:
    #             score = float(score)
    #         except:
    #             raise Exception('Found score value of {} in input'
    #                             .format(score))
    #
    #     if not tree:
    #         tree = s
    #
    #     if program == 'None':
    #         program = None
    #
    #     if taxon_set is not None:
    #         tree = cls(newick=tree, name=name, program=program, score=score,
    #                    taxon_set=taxon_set)
    #
    #     else:
    #         tree = cls(newick=tree, name=name, program=program, score=score)
    #
    #     return tree

    def write_to_file(
            self,
            outfile,
            metadata=False,
            scale=1,
            **kwargs
    ):
        """ Writes the tree to file. If metadata==True it writes the tree's
        likelihood score, name and generating program to the file inside a
        comment.
        Scale allows branch lengths to be scaled by a float, or suppressed if
        the scale value is zero. Negative values also suppress edge lengths.
        KWargs allowed are all the KWargs allowed by the dendropy tree.as_string
        method, including:
        suppress_edge_lengths, internal_labels, suppress_rooting,
            edge_label_compose_func...
        suppress_rooting is set to True by default
        """
        if scale <= 0:
            kwargs.update({'suppress_edge_lengths': True})
        else:
            l = lambda x: str(x.length * scale)
            kwargs.update({'edge_label_compose_func': l})
        if not 'suppress_rooting' in kwargs:
            kwargs.update({'suppress_rooting': True})

        with open(outfile, 'w') as writer:
            if metadata:
                writer.write(str(self))
            else:
                writeable = self.as_string('newick', **kwargs)
                writer.write(writeable + '\n')
        return outfile

    def __name_things(self):
        """ Easy names for debugging """
        edges = {}
        nodes = {None: 'root'}
        for n in self.postorder_node_iter():
            nodes[n] = '.'.join([str(x.taxon) for x in n.leaf_nodes()])
        for e in self.preorder_edge_iter():
            edges[e] = ' ---> '.join([nodes[e.tail_node], nodes[e.head_node]])

        r_edges = {value: key for key, value in edges.items()}
        r_nodes = {value: key for key, value in nodes.items()}
        return edges, nodes, r_edges, r_nodes

    # @classmethod
    # def new_tree_from_phyml_results(
    #         cls,
    #         tree_file,
    #         stats_file,
    #         name=None,
    #         program='phyml',
    # ):
    #     """ Given the usual phyml output files - xxx_phyml_tree.txt and
    #     xxx_phyml_stats.txt, instantiates a Tree from the information
    #     in the phyml files.
    #     TODO: refactor into phymlIO module """
    #     # newick score output program name
    #
    #     exit_ = False
    #     for f in (tree_file, stats_file):
    #         try:
    #             filecheck(f)
    #         except FileError, e:
    #             print(e)
    #             exit_ = True
    #
    #     if exit_:
    #         print('Results were not loaded')
    #         raise FileError()
    #
    #     if not name:
    #         name = cls.name_regex.search(tree_file).group(1)
    #     newick = open(tree_file).read()
    #     stats = open(stats_file).read()
    #     score = cls.score_regex.search(stats).group(0)
    #     score = float(score) if score else None
    #
    #     return cls(newick=newick, score=score, output=stats, program=program,
    #                name=name)name

    @classmethod
    def new_iterative_rtree(cls, nspecies):
        return RandomTree.new(nspecies)

    @classmethod
    def new_rtree(cls, nspecies=16, zero_root_height=True, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        tree = tg.rtree()
        if zero_root_height:
            tree.seed_node.edge_length = 0.0
        return tree

    @classmethod
    def new_coal(cls, nspecies=16, zero_root_height=True, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        tree = tg.coal()
        if zero_root_height:
            tree.seed_node.edge_length = 0.0
        return tree

    @classmethod
    def new_yule(cls, nspecies=16, zero_root_height=True, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        tree = tg.yule()
        if zero_root_height:
            tree.seed_node.edge_length = 0.0
        return tree

    def sample_gene_tree(self, **kwargs):
        tg = TreeGen(template=self)
        return tg.gene_tree(**kwargs)['gene_tree']


class RandomTree(object):
    def __init__(self, names=None):
        if names is None:
            self.label_generator = itertools.chain(_infinite_labels_generator(['l'], start=1))
            next(self.label_generator)
        else:
            self.label_generator = itertools.chain(_infinite_labels_generator(names, start=2))

        self.tree = Tree('({}:1,{}:1,{}:1):0'.format(self.next_label(), self.next_label(), self.next_label()))

    def next_label(self):
        return next(self.label_generator)

    def new_taxon_object(self):
        lab = self.next_label()
        tax = dendropy.Taxon(label=lab)
        return tax

    def add(self, edge):
        tail = edge.tail_node
        head = edge.head_node
        tail.remove_child(head)
        new_taxon = self.new_taxon_object()
        new_inner = tail.new_child(edge_length=1.0)
        new_inner.new_child(taxon=new_taxon, edge_length=1.0)
        new_inner.add_child(head, edge_length=1.0)

    def select(self):
        e, _ = self.tree.map_event_onto_tree()
        return e

    @classmethod
    def new(cls, n, names=None):
        rt = cls(names)
        for _ in range(n - 3):
            e = rt.select()
            rt.add(e)
        return rt.tree


class TreeGen(object):
    def __init__(
            self,
            nspecies=16,
            names=None,
            template=None,
            cf=False,
    ):

        """ Generates a new Tree using a coalescent process (coal method), a
        Yule pure-birth process (yule method), a random tree (rtree), or by
        sampling a gene tree from a template species tree using a constrained
        Kingman coalescent.

        nspecies = number of taxa in the tree
        names = a list of leaf names (names will be generated if not supplied)
        template = a template species tree for drawing gene trees
        cf = set to true to generate leaf names from the list of character names
        from Cannon Fodder """

        self.nspecies = nspecies
        if names is not None:
            g = _infinite_labels_generator(names, shuffle=False)
            self.names = list(itertools.islice(g, nspecies))
        if cf:
            g = _infinite_labels_generator(cfnames)
            self.names = list(itertools.islice(g, nspecies))
        else:
            g = itertools.chain(_infinite_labels_generator(['Sp'], start=1))
            next(g)
            self.names = list(itertools.islice(g, nspecies))

        if template and not isinstance(template, Tree):
            raise TypeError('template should be \'Tree\' object. Got',
                            type(template))
        self.template = template

    def coal(self):
        taxon_set = dendropy.TaxonSet(self.names)
        return cast(dendropy.treesim.pure_kingman(taxon_set))

    def gene_tree(
            self,
            scale_to=None,
            population_size=1,
            trim_names=True,
    ):
        """ Using the current tree object as a species tree, generate a gene
        tree using the constrained Kingman coalescent process from dendropy. The
        species tree should probably be a valid, ultrametric tree, generated by
        some pure birth, birth-death or coalescent process, but no checks are
        made. Optional kwargs are: -- scale_to, which is a floating point value
        to scale the total tree tip-to-root length to, -- population_size, which
        is a floating point value which all branch lengths will be divided by to
        convert them to coalescent units, and -- trim_names, boolean, defaults
        to true, trims off the number which dendropy appends to the sequence
        name """

        tree = self.template or self.yule()

        for leaf in tree.leaf_iter():
            leaf.num_genes = 1

        dfr = tree.seed_node.distance_from_root()
        dft = tree.seed_node.distance_from_tip()
        tree_height = dfr + dft

        if scale_to:
            population_size = tree_height / scale_to

        for edge in tree.preorder_edge_iter():
            edge.pop_size = population_size

        gene_tree = dendropy.treesim.constrained_kingman(tree)[0]

        if trim_names:
            for leaf in gene_tree.leaf_iter():
                leaf.taxon.label = leaf.taxon.label.replace('\'', '').split('_')[0]

        return {'gene_tree': tree.__class__(gene_tree.as_newick_string().strip(';') + ';'),
                'species_tree': tree}

    def rtree(self):
        m = self.yule()
        m.randomise_labels()
        return m.randomise_branch_lengths()

    def yule(self):
        taxon_set = dendropy.TaxonSet(self.names)
        return cast(dendropy.treesim.uniform_pure_birth(taxon_set))


cfnames = [
    'Jools', 'Jops', 'Stoo', 'Rj', 'Ubik', 'Cj', 'Chris', 'Pete',
    'Tadger', 'Hector', 'Elroy', 'Softy', 'Mac', 'Bomber', 'Stan', 'Tosh',
    'Brains', 'Norm', 'Buster', 'Spike', 'Browny', 'Murphy', 'Killer', 'Abdul',
    'Spotty', 'Goofy', 'Donald', 'Windy', 'Nifta', 'Denzil', 'Cedric', 'Alf',
    'Marty', 'Cecil', 'Wally', 'Pervy', 'Jason', 'Roy', 'Peewee', 'Arnie',
    'Lofty', 'Tubby', 'Porky', 'Norris', 'Bugsy', 'Greg', 'Gus', 'Ginger',
    'Eddy', 'Steve', 'Hugo', 'Zippy', 'Sonny', 'Willy', 'Mario', 'Luigi',
    'Bo', 'Johan', 'Colin', 'Queeny', 'Morgan', 'Reg', 'Peter', 'Brett',
    'Matt', 'Vic', 'Hut', 'Bud', 'Brad', 'Ashley', 'Les', 'Rex',
    'Louis', 'Pedro', 'Marco', 'Leon', 'Ali', 'Tyson', 'Tiger', 'Frank',
    'Reuben', 'Leyton', 'Josh', 'Judas', 'Aj', 'Lex', 'Butch', 'Bison',
    'Gary', 'Luther', 'Kermit', 'Brian', 'Ray', 'Freak', 'Leroy', 'Lee',
    'Banjo', 'Beaker', 'Basil', 'Bonzo', 'Kelvin', 'Ronnie', 'Rupert', 'Roo',
    'Dan', 'Jimmy', 'Bob', 'Don', 'Tommy', 'Eddie', 'Ozzy', 'Paddy',
    'Arnold', 'Tony', 'Teddy', 'Dom', 'Theo', 'Martin', 'Chunky', 'Jon',
    'Ben', 'Girly', 'Julian', 'Pizza', 'Ciaran', 'Jock', 'Gravy', 'Trendy',
    'Neil', 'Derek', 'Ed', 'Biff', 'Paul', 'Stuart', 'Randy', 'Loreta',
    'Suzie', 'Pumpy', 'Urmer', 'Roger', 'Pussy', 'Meat', 'Beefy', 'Harry',
    'Tiny', 'Howard', 'Morris', 'Thor', 'Rev', 'Duke', 'Micky', 'Chas',
    'Melony', 'Craig', 'Sidney', 'Parson', 'Rowan', 'Smelly', 'Dok', 'Stew',
    'Adrian', 'Pat', 'Iceman', 'Goose', 'Dippy', 'Viv', 'Fags', 'Bunty',
    'Noel', 'Bono', 'Edge', 'Robbie', 'Sean', 'Miles', 'Jimi', 'Gordon',
    'Val', 'Hobo', 'Fungus', 'Toilet', 'Lampy', 'Marcus', 'Pele', 'Hubert',
    'James', 'Tim', 'Saul', 'Andy', 'Silky', 'Simon', 'Handy', 'Sid',
    'George', 'Joff', 'Barry', 'Dick', 'Gil', 'Nick', 'Ted', 'Phil',
    'Woody', 'Wynn', 'Alan', 'Pip', 'Mickey', 'Justin', 'Karl', 'Maddog',
    'Horace', 'Harold', 'Gazza', 'Spiv', 'Foxy', 'Ned', 'Bazil', 'Oliver',
    'Rett', 'Scot', 'Darren', 'Moses', 'Noah', 'Seth', 'Buddah', 'Mary',
    'Pilot', 'Mcbeth', 'Mcduff', 'Belly', 'Mathew', 'Mark', 'Luke', 'John',
    'Aslam', 'Ham', 'Shem', 'Joshua', 'Jacob', 'Esaw', 'Omar', 'Enoch',
    'Obadia', 'Daniel', 'Samuel', 'Robbo', 'Joebed', 'Ismael', 'Isreal', 'Isabel',
    'Isarat', 'Monk', 'Blip', 'Bacon', 'Danube', 'Friend', 'Darryl', 'Izzy',
    'Crosby', 'Stills', 'Nash', 'Young', 'Cheese', 'Salami', 'Prawn', 'Radish',
    'Egbert', 'Edwy', 'Edgar', 'Edwin', 'Edred', 'Eggpie', 'Bros', 'Sonic',
    'Ziggy', 'Alfred', 'Siggy', 'Hilda', 'Snell', 'Sparks', 'Spook', 'Topcat',
    'Benny', 'Dibble', 'Benker', 'Dosey', 'Beaky', 'Joist', 'Pivot', 'Tree',
    'Bush', 'Grass', 'Seedy', 'Tin', 'Rollo', 'Zippo', 'Nancy', 'Larry',
    'Iggy', 'Nigel', 'Jamie', 'Jesse', 'Leo', 'Virgo', 'Garth', 'Fidel',
    'Idi', 'Che', 'Kirk', 'Spock', 'Maccoy', 'Chekov', 'Uhura', 'Bones',
    'Vulcan', 'Fester', 'Jethro', 'Jimbob', 'Declan', 'Dalek', 'Hickey', 'Chocco',
    'Goch', 'Pablo', 'Renoir', 'Rolf', 'Dali', 'Monet', 'Manet', 'Gaugin',
    'Chagal', 'Kid', 'Hully', 'Robert', 'Piers', 'Raith', 'Jeeves', 'Paster',
    'Adolf', 'Deiter', 'Deni', 'Zark', 'Wizkid', 'Wizard', 'Iain', 'Kitten',
    'Gonner', 'Waster', 'Loser', 'Fodder',
]
