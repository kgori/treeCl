#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import next
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div

# standard library
import itertools
import random

# third party
import dendropy as dpy
import numpy as np
from tree_distance import PhyloTree

# treeCl
from .errors import optioncheck
from .constants import ISPY3
from .utils import fileIO, weighted_choice
from .utils.decorators import lazyprop
from .utils.math import truncated_exponential

import logging
logger = logging.getLogger(__name__)


def cast(dendropy_tree):
    """ Cast dendropy.Tree instance as Tree instance """
    return Tree(dendropy_tree.as_string('newick', suppress_rooting=True) + ';')


def _infinite_labels_generator(labels, start=2, shuffle=True):
    l = len(labels)
    loop1 = random.sample(labels, l) if shuffle else labels
    return itertools.chain.from_iterable([loop1, ('{}{}'.format(x, y) for x, y in
                                                  zip(itertools.cycle(labels),
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
        root_child_nodes = self.tree._tree.seed_node.child_nodes()
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
        self.tree._tree.prune_subtree(n, suppress_unifurcations=False)
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
        new.add_child(h)
        h.edge.length=length
        new.add_child(node)
        self.tree._dirty = True
        self.tree._tree.encode_bipartitions(suppress_unifurcations=True)

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

        starting_length = self.tree._tree.length()

        excl = [self.tree._tree.seed_node.edge]  # exclude r
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
            if set(self.tree._tree.preorder_edge_iter()) - set(excl) == set([]):
                print(repr(self.tree))
                print(self.tree._tree.as_ascii_plot())
                # print(edges[prune_edge])
                raise Exception('No non-sibling sprs available')

        regraft_edge, l2 = self.tree.map_event_onto_tree(excl)

        # edges, nodes, redges, rnodes = self.tree._name_things()
        # print(edges[prune_edge], l1, edges[regraft_edge], l2)
        self.spr(prune_edge, l1, regraft_edge, l2)
        if rescale:
            self.tree.scale(starting_length / self.tree.length())
            self.tree._dirty = True


class LGT(object):
    def __init__(self, tree):
        self.SPR = SPR(tree)
        self.tree = self.SPR.tree
        try:
            self.tree._tree.calc_node_ages()
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

        matching_edges = self.tree._tree.preorder_edge_iter(edge_matches_time)
        return list(matching_edges)

    def rlgt(self, time=None, disallow_sibling_lgts=False):
        self.tree._tree.calc_node_ages()
        excl = [self.tree._tree.seed_node.edge]

        if time is None:
            if disallow_sibling_lgts:
                self.add_single_node()
                children = self.tree._tree.seed_node.child_nodes()
                excl.extend([n.edge for n in children])
                time = self.get_time(excl)
                print('time = {0}'.format(time))
                self.tree._tree.encode_bipartitions()
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
        self.tree._tree.calc_node_ages()

    def add_single_node(self):
        cn = self.tree._tree.seed_node.child_nodes()
        el = lambda n: n.edge_length
        sh = min(cn, key=el)
        lo = max(cn, key=el)
        new = self.tree._tree.seed_node.new_child(edge_length=sh.edge_length)
        self.tree._tree.prune_subtree(lo, suppress_unifurcations=False)
        lo.edge_length -= sh.edge_length
        new.add_child(lo)
        self.tree._tree.encode_bipartitions(suppress_unifurcations=False)
        self.tree._tree.calc_node_ages()
        self.tree._dirty = True


class NNI(object):
    def __init__(self, tree):
        self.tree = tree

    def _validate(self):
        excludes = [self.tree._tree.seed_node] + self.tree._tree.leaf_nodes()
        if self.tree.rooted:
            child_a, child_b = self.tree._tree.seed_node.child_nodes()
        if child_a in excludes:
            excludes.append(child_b)
        if child_b in excludes:
            excludes.append(child_a)
        self.valid_nodes = set([n for n in self.tree._tree.nodes() if not n in excludes])

    def choose_node(self, use_weighted_choice=False, transform=None):
        self._validate()
        if use_weighted_choice:
            weights = np.array([n.edge.length for n in self.valid_nodes])
            if any(weight is None for weight in weights):
                logger.debug('Not all weights were valid: {}'.format(weights))
                weights = np.array([1.0 for n in self.valid_nodes])
            logger.debug('Weights (weighted choice=True): {}'.format(weights))
            if transform is not None:
                weights = transform(weights)
                logger.debug('Weights (transform=not None): {}'.format(weights))
        else:
            weights = np.array([1.0 for n in self.valid_nodes])
            logger.debug('Weights (weighted choice=False): {}'.format(weights))

        return weighted_choice(list(zip(self.valid_nodes, weights)))

    def get_exchangeable_nodes(self, n):
        """
            A      C    | Subtrees A, B, C and D are the exchangeable nodes
             \    /     | around the edge headed by n
              -->n      | The NNI exchanges either A or B with either C or D
             /    \
            B      D

            A      C                        C      A    | Subtree A is exchanged
             \    /        +NNI(A,C)         \    /     | with subtree C.
              -->n        ==========>         -->n
             /    \                          /    \
            B      D                        B      D
        """
        parent = n.parent_node
        a, b = random.sample(n.child_nodes(), 2)
        if parent.parent_node is None:
            if self.tree.rooted:
                c, d = random.sample(n.sister_nodes()[0].child_nodes(), 2)
            else:
                c, d = random.sample(n.sister_nodes(), 2)
        else:
            c = random.choice(n.sister_nodes())
            d = random.choice(parent.sister_nodes())

        return a, b, c, d

    def do_nni(self, node_1, node_2):
        parent_1 = node_1.parent_node
        parent_2 = node_2.parent_node
        parent_1.remove_child(node_1)
        parent_2.remove_child(node_2)
        parent_1.add_child(node_2)
        parent_2.add_child(node_1)
        self.tree._tree.encode_bipartitions()

    def rnni(self, use_weighted_choice=False, transform=None):
        n = self.choose_node(use_weighted_choice, transform)
        a, b, c, d = self.get_exchangeable_nodes(n)
        self.do_nni(random.choice([a, b]), random.choice([c, d]))


def collapse(t, threshold=None, keep_lengths=True, support_key=None, length_threshold=0.0):

    to_collapse = []
    for node in t._tree.postorder_node_iter():
        if node.is_leaf():
            if node.edge_length < length_threshold:
                node.edge_length = 0
            continue
        if node is t.seed_node:
            continue
        if threshold is not None:
            try:
                if support_key:
                    support = float(node.annotations.get_value(support_key))
                    node.label = support
                else:
                    support = float(node.label)
            except TypeError as e:
                raise SupportValueError('Inner node with length {} has no support value'.format(node.edge_length), e)
            except ValueError as e:
                raise SupportValueError(
                    'Inner node with length {} has a non-numeric support value {}'.format(node.edge_length), e)
            if support < threshold:
                to_collapse.append(node.edge)
        if node.edge_length < length_threshold:
            to_collapse.append(node.edge)

    for edge in to_collapse:
        if keep_lengths:
            for child in edge.head_node.child_nodes():
                child.edge.length += edge.length
        edge.collapse()
    return t


class UltrametricNNI(NNI):
    def __init__(self, tree):
        super(UltrametricNNI, self).__init__(tree)

    def _make_ultrametric(self):
        leaves = list(self.tree._tree.leaf_node_iter())
        root_tip_dists = [leaf.distance_from_root() for leaf in leaves]
        mean_tip_dist = np.mean(root_tip_dists)
        for dist, leaf in zip(root_tip_dists, leaves):
            leaf.edge.length += (mean_tip_dist - dist)

    def _validate(self):
        super(UltrametricNNI, self)._validate()
        self._make_ultrametric()
        self.tree._tree.calc_node_ages()

    def do_nni(self, node1, node2, node3, node4):
        pass


class ILS(object):
    def __init__(self, tree):
        self.minlen=0
        self.tree = tree
        self._validate()

    def _make_ultrametric(self):
        leaves = list(self.tree._tree.leaf_node_iter())
        root_tip_dists = [leaf.distance_from_root() for leaf in leaves]
        mean_tip_dist = np.mean(root_tip_dists)
        for dist, leaf in zip(root_tip_dists, leaves):
            leaf.edge.length += (mean_tip_dist - dist)

    def _break_ties(self):
        collapse(self.tree, keep_lengths=True, length_threshold=self.minlen)
        self.tree._tree.resolve_polytomies()

    def _validate(self):
        for edge in self.tree._tree.preorder_edge_iter():
            if np.isnan(edge.length) or edge.length <= self.minlen:
                edge.length = self.minlen
        self._break_ties()
        self._make_ultrametric()
        self.tree._tree.calc_node_ages()
        excludes = [self.tree._tree.seed_node] + self.tree._tree.seed_node.child_nodes() + self.tree._tree.leaf_nodes()
        self.valid_nodes = self.tree._tree.nodes(filter_fn=lambda x: not x in excludes)

    def choose_node(self, use_weighted_choice=False, transform=None):
        self._validate()
        if use_weighted_choice:
            weights = np.array([n.edge.length for n in self.valid_nodes])
            if any(weight is None for weight in weights):
                logger.debug('Not all weights were valid: {}'.format(weights))
                weights = np.array([1.0 for n in self.valid_nodes])
            logger.debug('Weights (weighted choice=True): {}'.format(weights))
            if transform is not None:
                weights = transform(weights)
                logger.debug('Weights (transform=not None): {}'.format(weights))
        else:
            weights = np.array([1.0 for n in self.valid_nodes])
            logger.debug('Weights (weighted choice=False): {}'.format(weights))

        return weighted_choice(list(zip(self.valid_nodes, weights)))

    def get_matching_edge(self, starting_node, time):
        def edge_matches_time(edge):
            if edge.tail_node is None:
                return False
            return edge.head_node.age < time < edge.tail_node.age

        if time > starting_node.parent_node.age:
            for node in starting_node.parent_node.ancestor_iter():
                if edge_matches_time(node.edge):
                    return node.edge
        else:
            sister = starting_node.sister_nodes()[0].edge
            if edge_matches_time(sister):
                return sister
            else:
                raise ValueError('No matching edge was found')


    def ils(self, node, sorting_times=None, force_topology_change=True):
        """
        A constrained and approximation of ILS using nearest-neighbour interchange

        Process
        -------
        A node with at least three descendents is selected from an ultrametric tree
        (node '2', below)

                   ---0--...              ---0--...           ---0--...
                  |        |             |        |        --1--      |
                  |        R           --1--      R       |     |     R
       age        |                   |     |            -2-    |
        ^         |                   |     |           |   |   |
        |       --1--                -2-    |           |   |   |
        |      |     |       or     |   |   |     or    |   |   |
        |      |     |              |   |   |           |   |   |
        |     -2-    |              |   |   |           |   |   |
        |    |   |   |              |   |   |           |   |   |
        |    A   B   C              C   B   A           A   C   B

        Nodes 'A', 'B' and 'C' are rearranged into one of the three configurations
        [(A, B), C], [A, (B, C)], [(A, C), B]

        Nodes 1 and 2 are slid further up the tree, but no further than node 0
        (this is why it's a constrained version), by an amount drawn from a
        truncated exponential distribution.

        This is approximately corresponds to the case where A and B failed to
        coalesce in the branch 1->2, so they coalesce with C in the branch
        0 -> 1 instead
        """
        # node = '2', par = '1', gpar = '0' -- in above diagram
        n_2 = node
        n_1 = n_2.parent_node
        if n_1 == self.tree._tree.seed_node:
            logger.warn('Node 1 is the root - calling again on child')
            self.ils(n_2.child_nodes())
        n_0 = n_1.parent_node
        a, b = node.child_nodes()
        c, = node.sister_nodes()

        ages = [a.age, b.age, c.age, n_2.age, n_1.age, n_0.age]

        # Do topology changes
        if force_topology_change:
            swap_mode = random.choice([1, 2])
        else:
            swap_mode = random.choice([0, 1, 2])

        if swap_mode == 1:
            # Exchange 'a' and 'c'
            n_2.remove_child(a)
            n_1.remove_child(c)
            n_2.add_child(c)
            n_1.add_child(a)

        elif swap_mode == 2:
            # Exchange 'b' and 'c'
            n_2.remove_child(b)
            n_1.remove_child(c)
            n_2.add_child(c)
            n_1.add_child(b)

        # Do branch length adjustments
        # Bounds - between node 0 (upper) and node 1 (lower)
        min_unsorted_age = n_1.age
        max_unsorted_age = n_0.age
        if sorting_times is None:
            sorting_times = truncated_exponential(max_unsorted_age-min_unsorted_age,
                                              scale=0.1*(max_unsorted_age-min_unsorted_age),
                                              sample_size=2) # E(t) = n(n-1)/2, n = 3
            sorting_times += min_unsorted_age
            sorting_times = np.array([min_unsorted_age, ages[3]])

        # Adjust node 1 edge length
        new_n1_age = max(sorting_times)
        prev_age = ages[4]
        slide = (new_n1_age - prev_age)
        if slide < 1e-6:
            slide = 0
            new_n1_age = prev_age
        n_1.edge.length -= slide
        n_2.edge.length += slide

        # Adjust node 2 edge length
        new_n2_age = min(sorting_times)
        prev_age = ages[3]
        slide = (new_n2_age - prev_age)
        if slide < 1e-6:
            slide = 0
            new_n2_age = prev_age
        n_2.edge.length -= slide

        # Adjust a, b and c edge lengths
        if swap_mode == 0:
            a.edge.length = (new_n2_age - ages[0])
            b.edge.length = (new_n2_age - ages[1])
            c.edge.length = (new_n1_age - ages[2])

        elif swap_mode == 1:
            a.edge.length = (new_n1_age - ages[0])
            b.edge.length = (new_n2_age - ages[1])
            c.edge.length = (new_n2_age - ages[2])

        else:
            a.edge.length = (new_n2_age - ages[0])
            b.edge.length = (new_n1_age - ages[1])
            c.edge.length = (new_n2_age - ages[2])

        # used to be .reindex_taxa() before dendropy 4.
        # migrate_taxon_namespace is recommended migrated function,
        # but not sure if its even needed anymore.
        self.tree._tree.migrate_taxon_namespace(self.tree._tree.taxon_namespace)

        self.tree._tree.encode_bipartitions()
        self._validate()
        logger.debug(self.tree)

    def rils(self, use_weighted_choice=True, transform=None):
        n = self.choose_node(use_weighted_choice, transform)
        logger.debug('Chosen node = {} age = {} parent age = {}'.format([leaf.taxon.label for leaf in n.leaf_nodes()], n.age, n.parent_node.age))
        self.ils(n)

    # def ils_(self, node, sorting_times=None):
    #     unsorted_descendants = node.child_nodes()
    #     logger.info('Child 1 = {} age = {}'.format([leaf.taxon.label for leaf in unsorted_descendants[0].leaf_nodes()], unsorted_descendants[0].age))
    #     logger.info('Child 1 = {} age = {}'.format([leaf.taxon.label for leaf in unsorted_descendants[1].leaf_nodes()], unsorted_descendants[1].age))

    #     min_unsorted_age = max(node.age, node.sister_nodes()[0].age)
    #     logger.debug('Node age = {}, sister age = {}'.format(node.age, node.sister_nodes()[0].age))
    #     max_unsorted_age = self.tree.seed_node.age

    #     if sorting_times is None:
    #         sorting_times = truncated_exponential(max_unsorted_age-min_unsorted_age,
    #                                           scale=0.5*(max_unsorted_age-min_unsorted_age),
    #                                           sample_size=2) # E(t) = n(n-1)/2, n = 2
    #         sorting_times += min_unsorted_age

    #     if np.any(sorting_times > max_unsorted_age): logger.error('Sorting times too large: {}'.format(sorting_times))
    #     logger.info('Min/Max ages = {} {}'.format(min_unsorted_age, max_unsorted_age))
    #     logger.info('Sorting occurs = {} {}'.format(*sorting_times))

    #     random.shuffle(unsorted_descendants)
    #     c1, c2 = unsorted_descendants
    #     time1 = max(sorting_times)
    #     time2 = min(sorting_times)
    #     donor1 = self.get_matching_edge(node, time1)
    #     donor2 = self.get_matching_edge(node, time2)

    #     logger.info('Stage 0 - initial tree')
    #     logger.info('Stage 1 - remove c1')
    #     self.tree.print_plot(plot_metric='length')
    #     node.remove_child(c1)

    #     logger.info('Stage 2 - remove c2')
    #     self.tree.print_plot(plot_metric='length')
    #     node.remove_child(c2)

    #     c1.edge.length = time1 - c1.age
    #     c2.edge.length = time2 - c2.age

    #     logger.info('Stage 3 - regraft c1 at time={}'.format(time1))
    #     self.tree.print_plot(plot_metric='length')
    #     self.SPR.regraft(donor1, c1, time1 - donor1.head_node.age)

    #     logger.info('Stage 4 - regraft c2 at time={}'.format(time2))
    #     self.tree.print_plot(plot_metric='length')
    #     self.SPR.regraft(donor2, c2, time2 - donor2.head_node.age)

    #     self.tree.print_plot(plot_metric='length')
    #     self.tree.prune_subtree(node, suppress_unifurcations=True)

    #     logger.info('Stage Final')
    #     self.tree.print_plot(plot_metric='length')
    #     self.tree.encode_bipartitions()
    #     self.tree.reindex_taxa()
    #     self.tree.calc_node_ages()
    #     self._validate()


class NNI2(object):
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
        if not self.tree._tree.seed_node == t:
            original_seed = self.tree._tree.seed_node
            self.tree._tree.reseed_at(t)
        else:
            original_seed = None
        head_children = h.child_nodes()
        tail_children = list(set(t.child_nodes()) - {h})  # See N1
        if original_seed:
            self.tree._tree.reseed_at(original_seed)

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

        original_seed = self.tree._tree.seed_node
        head = edge.head_node
        tail = edge.tail_node
        self.tree._tree.reseed_at(tail)
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
        self.tree._tree.reseed_at(original_seed)
        self.tree._tree.encode_bipartitions()
        self.tree._dirty = True

    def reroot_tree(self):
        if self.reroot and self.rooting_info is not None:
            self.tree._tree.reroot_at_edge(*self.rooting_info)
            self.tree._tree.encode_bipartitions()
            self.tree._dirty = True
        return self.tree

    def rnni(self, use_weighted_choice=False, invert_weights=False):
        """
        Apply a random NNI operation at a randomly selected edge
        The edge can be chosen uniformly, or weighted by length --
        invert_weights favours short edges.
        """
        if use_weighted_choice:
            leaves = list(self.tree._tree.leaf_edge_iter())
            e, _ = self.tree.map_event_onto_tree(excluded_edges=leaves, invert_weights=invert_weights)
        else:
            e = random.choice(self.tree.get_inner_edges())
        children = self.get_children(e)
        h = random.choice(children['head'])
        t = random.choice(children['tail'])
        self.nni(e, h, t)


class Tree(object):
    """ Tree data structure, wraps dendropy Tree class
    """

    def __init__(
            self,
            newick=None,
            name=None,
            **kwargs
    ):

        if newick:
            self._tree = dpy.Tree.get_from_string(newick, 'newick', preserve_underscores=True, **kwargs)
            if self.rooted:
                self._tree.is_rooted = True
                self._tree.encode_bipartitions()
        else:
            self._tree = dpy.Tree(**kwargs)

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
        return len(self._tree.leaf_nodes())

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
        return set([n.taxon.label for n in self._tree.leaf_nodes()])

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
        suppress_internal_node_labels [True/False]
            - turn on/off bootstrap labels
        suppress_rooting [True/False]
            - turn on/off [&U] or [&R] rooting
              state labels
        edge_label_compose_func
            - function to convert edge lengths:
              takes edge as arg, returns string
        """
        n = self._tree.as_string('newick',
                                 suppress_rooting=True,
                                 suppress_internal_node_labels=True)
        if n:
            return n.strip(';\n') + ';'
        return n

    @property
    def phylotree(self):
        """
        Get the c++ PhyloTree object corresponding to this tree.
        :return: PhyloTree instance
        """
        if not self._phylotree or self._dirty:
            try:
                if ISPY3:
                    self._phylotree = PhyloTree(self.newick.encode(), self.rooted)
                else:
                    self._phylotree = PhyloTree(self.newick, self.rooted)
            except ValueError:
                logger.error('Couldn\'t convert to C++ PhyloTree -- are there bootstrap values?')
            self._dirty = False
        return self._phylotree

    @property
    def seed_node(self):
        return self._tree.seed_node

    @newick.setter
    def newick(self, newick_string):
        if self.newick:
            print('Newick string already loaded: {0}'.format(self.newick))
            return
        self._tree = dpy.Tree.get_from_string(newick_string, 'newick')

    @property
    def rooted(self):
        """ Predicate testing for rootedness by checking for a bifurcation
        at the root. """
        return len(self._tree.seed_node.child_nodes()) == 2 if self.newick else None

    @classmethod
    def bifurcate_base(cls, newick):
        """ Rewrites a newick string so that the base is a bifurcation
        (rooted tree) """
        t = cls(newick)
        t._tree.resolve_polytomies()
        return t.newick

    @classmethod
    def trifurcate_base(cls, newick):
        """ Rewrites a newick string so that the base is a trifurcation
        (usually means an unrooted tree) """
        t = cls(newick)
        t._tree.deroot()
        return t.newick

    def copy(self):
        """ Returns an independent copy of self """
        return self.__class__(self.newick)

    def deroot(self):
        """ Unroot the tree, inplace """
        self._tree.deroot()

    def get_inner_edges(self):
        """ Returns a list of the internal edges of the tree. """
        inner_edges = [e for e in self._tree.preorder_edge_iter() if e.is_internal()
                       and e.head_node and e.tail_node]
        return inner_edges

    def get_nonroot_edges(self):
        return [e for e in self._tree.preorder_edge_iter()
                if e.head_node and e.tail_node]

    def intersection(self, other):
        """ Returns the intersection of the taxon sets of two Trees """
        taxa1 = self.labels
        taxa2 = other.labels
        return taxa1 & taxa2

    def map_event_onto_tree(self, excluded_edges=None, invert_weights=False):
        edge_list = list(self._tree.preorder_edge_iter())
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
        if invert_weights:
            lengths = 1/lengths
        cumulative_lengths = lengths.cumsum()
        rnum = np.random.random() * cumulative_lengths[-1]
        index = cumulative_lengths.searchsorted(rnum)
        chosen_edge = edge_list[index]
        from_head_length = cumulative_lengths[index] - rnum
        return chosen_edge, from_head_length

    def multifurcate(self, threshold=1e-06, update_splits=True):
        for edge in self._tree.postorder_edge_iter():
            if edge.is_internal():
                if edge.length <= threshold:
                    edge.collapse()
                    self._dirty = True
        if update_splits:
            self._tree.encode_bipartitions()

    def ntaxa(self):
        return len(self)

    def pairdist(self, taxon_label1, taxon_label2):
        if self.patristic is None:
            print('Error calculating patristic distances - maybe this '
                  'tree has no branch lengths?')
            return

        leaf1 = self._tree.find_node_with_taxon_label(taxon_label1)
        leaf2 = self._tree.find_node_with_taxon_label(taxon_label2)

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
            pdm = dpy.calculate.treemeasure.PatristicDistanceMatrix(self._tree)
        except TypeError:
            pdm = None
        return pdm

    def postorder(self, skip_seed=False):
        """
        Return a generator that yields the nodes of the tree in postorder.
        If skip_seed=True then the root node is not included.
        """
        for node in self._tree.postorder_node_iter():
            if skip_seed and node is self._tree.seed_node:
                continue
            yield node

    def preorder(self, skip_seed=False):
        """
        Return a generator that yields the nodes of the tree in preorder.
        If skip_seed=True then the root node is not included.
        """
        for node in self._tree.preorder_node_iter():
            if skip_seed and node is self._tree.seed_node:
                continue
            yield node

    def prune_to_subset(self, subset, inplace=False):
        """ Prunes the Tree to just the taxon set given in `subset` """
        if not subset.issubset(self.labels):
            print('"subset" is not a subset')
            return
        if not inplace:
            t = self.copy()
        else:
            t = self
        t._tree.retain_taxa_with_labels(subset)
        t._tree.encode_bipartitions()
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

        for n in t._tree.preorder_node_iter():
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

        names = list(t.labels)
        random.shuffle(names)
        for l in t._tree.leaf_node_iter():
            l.taxon._label = names.pop()
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
        root_edge = self._tree.seed_node.edge
        lengths = dict([(edge, edge.length) for edge
                        in self._tree.seed_node.incident_edges() if edge is not root_edge])
        self._tree.deroot()
        reroot_edge = (set(self._tree.seed_node.incident_edges())
                       & set(lengths.keys())).pop()
        self._tree.encode_bipartitions()
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
            for node in self._tree.preorder_node_iter():
                node.rate = root_rate
            return

        for node in self._tree.preorder_node_iter():
            if node == self._tree.seed_node:
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

        for node in self._tree.preorder_node_iter():
            if node == self._tree.seed_node:
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

    def rnni(self, times=1, **kwargs):
        """ Applies a NNI operation on a randomly chosen edge.
        keyword args: use_weighted_choice (True/False) weight the random edge selection by edge length
                      transform (callable) transforms the edges using this function, prior to weighted selection
        """

        nni = NNI(self.copy())
        for _ in range(times):
            nni.rnni(**kwargs)
        # nni.reroot_tree()
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
        t._tree.scale_edges(factor)
        t._dirty = True
        return t

    def strip(self, inplace=False):
        """ Sets all edge lengths to None """
        if not inplace:
            t = self.copy()
        else:
            t = self
        for e in t._tree.preorder_edge_iter():
            e.length = None
        t._dirty = True
        return t

    def translate(self, dct):
        """
        Translate leaf names using a dictionary of names
        :param dct: Dictionary of current names -> updated names
        :return: Copy of tree with names changed
        """
        new_tree = self.copy()
        for leaf in new_tree._tree.leaf_node_iter():
            curr_name = leaf.taxon.label
            leaf.taxon.label = dct.get(curr_name, curr_name)
        return new_tree

    def _name_things(self):
        """ Easy names for debugging """
        edges = {}
        nodes = {None: 'root'}
        for n in self._tree.postorder_node_iter():
            nodes[n] = '.'.join([str(x.taxon) for x in n.leaf_nodes()])
        for e in self._tree.preorder_edge_iter():
            edges[e] = ' ---> '.join([nodes[e.tail_node], nodes[e.head_node]])

        r_edges = {value: key for key, value in edges.items()}
        r_nodes = {value: key for key, value in nodes.items()}
        return edges, nodes, r_edges, r_nodes

    @classmethod
    def new_iterative_rtree(cls, nspecies, **kwargs):
        return RandomTree.new(nspecies, **kwargs)

    @classmethod
    def new_rtree(cls, nspecies=16, zero_root_height=True, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        tree = tg.rtree()
        if zero_root_height:
            tree._tree.seed_node.edge_length = 0.0
        return tree

    @classmethod
    def new_coal(cls, nspecies=16, zero_root_height=True, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        tree = tg.coal()
        if zero_root_height:
            tree._tree.seed_node.edge_length = 0.0
        return tree

    @classmethod
    def new_yule(cls, nspecies=16, zero_root_height=True, **kwargs):
        tg = TreeGen(nspecies, **kwargs)
        tree = tg.yule()
        if zero_root_height:
            tree._tree.seed_node.edge_length = 0.0
        return tree

    def sample_gene_tree(self, **kwargs):
        tg = TreeGen(template=self)
        return tg.gene_tree(**kwargs)['gene_tree']


class RandomTree(object):
    def __init__(self, names=None, rooted=False):
        if names is None:
            self.label_generator = itertools.chain(_infinite_labels_generator(['l'], start=1))
            next(self.label_generator)
        else:
            self.label_generator = itertools.chain(_infinite_labels_generator(names, start=2))

        if rooted:
            self.tree = Tree('({}:1,{}:1):0;'.format(self.next_label(), self.next_label()))
        else:
            self.tree = Tree('({}:1,{}:1,{}:1):0;'.format(self.next_label(), self.next_label(), self.next_label()))

    def next_label(self):
        return next(self.label_generator)

    def new_taxon_object(self):
        lab = self.next_label()
        tax = dpy.Taxon(label=lab)
        return tax

    def add(self, edge):
        tail = edge.tail_node
        head = edge.head_node
        tail.remove_child(head)
        new_taxon = self.new_taxon_object()
        new_inner = tail.new_child(edge_length=1.0)
        new_inner.new_child(taxon=new_taxon, edge_length=1.0)
        new_inner.add_child(head)
        head.edge_length=1.0

    def select(self):
        e, _ = self.tree.map_event_onto_tree()
        return e

    @classmethod
    def new(cls, n, names=None, rooted=False):
        rt = cls(names, rooted)
        present = 2 if rooted else 3
        for _ in range(n - present):
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
        elif cf:
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
        taxon_set = dpy.TaxonNamespace(self.names)
        return cast(dpy.simulate.treesim.pure_kingman_tree(taxon_set))

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

        for leaf in tree._tree.leaf_node_iter():
            leaf.num_genes = 1

        dfr = tree._tree.seed_node.distance_from_root()
        dft = tree._tree.seed_node.distance_from_tip()
        tree_height = dfr + dft

        if scale_to:
            population_size = tree_height / scale_to

        for edge in tree._tree.preorder_edge_iter():
            edge.pop_size = population_size

        gene_tree = dpy.simulate.treesim.constrained_kingman_tree(tree._tree)[0]

        if trim_names:
            for leaf in gene_tree.leaf_node_iter():
                leaf.taxon.label = leaf.taxon.label.replace('\'', '').split('_')[0]

        # Dendropy changed its API
        return {'gene_tree': tree.__class__(gene_tree.as_string('newick', suppress_rooting=True).strip(';\n') + ';'),
                'species_tree': tree}


    def rtree(self):
        m = self.yule()
        m.randomise_labels()
        return m.randomise_branch_lengths()

    def yule(self):
        taxon_set = dpy.TaxonNamespace(self.names)
        return cast(dpy.simulate.treesim.uniform_pure_birth_tree(taxon_set))


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
