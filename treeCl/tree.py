#!/usr/bin/env python

import  random
import  re
import  dendropy        as      dpy
from    numpy.random    import  gamma
from    errors          import  FileError, filecheck_and_raise
import  utils.dpy
from    externals       import  GTP


class TreeManip(object):

    def __init__(self, tree):
        if not isinstance(tree, Tree):
            raise TypeError('TreeManip class should be initialised with \'Tree\' object. Got'
                            , type(tree))
        self.tree = tree.copy()


    def __str__(self):
        return 'TreeManip object with tree:\n' + str(self.tree)

    def convert_to_dendropy_tree(self):
        """Takes Treee object, returns dendropy.Tree object"""

        return utils.dpy.convert_to_dendropy_tree(self.tree)

    def dendropy_as_newick(self, dpy_tree):
        return utils.dpy.convert_dendropy_to_newick(dpy_tree)

    def randomise_branch_lengths(
        self,
        i=(1, 1),
        l=(1, 1),
        distribution_func=gamma,
        output_format=5,
        ):
        """ i and l are tuples describing the parameters of the distribution
        function for inner edges and leaves, respectively. distribution_func is
        a function generating samples from a probability distribution (eg gamma,
        normal ...) """

        t = self.convert_to_dendropy_tree()
        for n in t.preorder_node_iter():
            if n.is_internal():
                n.edge.length = max(0, distribution_func(*i))
            else:
                n.edge.length = max(0, distribution_func(*l))
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def randomise_labels(self):
        t = self.convert_to_dendropy_tree()
        names = [l.taxon.label for l in t.leaf_iter()]
        names_copy = names[:]
        random.shuffle(names_copy)
        for l in t.leaf_iter():
            l.taxon.label = names_copy.pop()
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def scale(self, scaling_factor):

        t = self.convert_to_dendropy_tree()
        for e in t.preorder_edge_iter():
            if e.length:
                e.length *= scaling_factor
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def strip(self):

        t = self.convert_to_dendropy_tree()
        for e in t.preorder_edge_iter():
            e.length = None
        tree_copy = self.tree.copy()
        tree_copy.newick = self.dendropy_as_newick(t)
        self.tree = tree_copy
        return self.tree

    def nni(self):
        """Function to perform a random nearest-neighbour interchange on a tree
        using Dendropy
        
        The dendropy representation of a tree is as if rooted (even when it's
        not) The node which acts like the root is called the seed node, and this
        can sit in the middle of an edge which would be a leaf edge in the
        unrooted tree. NNI moves are only eligible on internal edges, so we need
        to check if the seed node is sat on a real internal edge, or a fake
        one."""

        tree = self.convert_to_dendropy_tree()
        seed = tree.seed_node
        resolved = False                    
        if len(seed.child_nodes()) > 2:     # If root is multifurcating we
            print 'Resolve root trisomy'    # need to force it to bifurcate,
            tree.resolve_polytomies()       # and re-establish multifurcation
            resolved = true                 # later

        # Make a list of internal edges not including the root edge
        edge_list = list(tree.preorder_edge_iter(lambda edge: \
                         (True if edge.is_internal() and edge.head_node != seed
                         and edge.tail_node != seed else False)))

        # Test whether root edge is eligible (i.e., the edge is not a
        # leaf when the tree is unrooted). If the test passes, add 'root'
        # to the list of eligible edges
        if not any([x.is_leaf() for x in seed.child_nodes()]):
            edge_list += ['root']

        chosen_edge = random.choice(edge_list)  # Choose the edge around which
        print chosen_edge                       # to do the NNI

        # Reroot at chosen edge, if necessary
        if chosen_edge != 'root':
            tree.reroot_at_edge(chosen_edge, length1=chosen_edge.length / 2,
                                length2=chosen_edge.length / 2,
                                delete_outdegree_one=False)
            root = tree.seed_node
        else:
            root = seed

        # To do the swap: find the nodes on either side of root
        (child_left, child_right) = root.child_nodes()

        # Pick a child node from each of these
        neighbour1 = random.choice(child_left.child_nodes())
        neighbour2 = random.choice(child_right.child_nodes())

        # Prune the chosen nearest neighbours - but don't
        # do any tree structure updates
        tree.prune_subtree(neighbour1, update_splits=False,
                           delete_outdegree_one=False)
        tree.prune_subtree(neighbour2, update_splits=False,
                           delete_outdegree_one=False)

        # Reattach the pruned neighbours to the opposite side
        # of the tree
        child_left.add_child(neighbour2)
        child_right.add_child(neighbour1)

        # Reroot the tree using the original seed node, and
        # update splits
        if not chosen_edge == 'root':
            tree.reroot_at_node(seed, update_splits=True)
        else:
            tree.update_splits()

        if resolved:
            print 'Reinstating root trisomy'
            tree.deroot()

        newick = self.dendropy_as_newick(tree)
        if tree.is_rooted:
            newick = '[&R] ' + newick

        tree_copy = self.tree.copy() 
        tree_copy.newick = newick
        self.tree = tree_copy
        return self.tree

    def spr(self, time=None, disallow_sibling_SPRs=False, verbose=False):

        def _get_blocks(tree, include_leaf_nodes=True):
            dists = []
            blocks = {}
            if include_leaf_nodes:
                iterator = tree.preorder_node_iter()
            else:
                iterator = tree.preorder_internal_node_iter()
            for n in iterator:
                node_height = n.distance_from_root()
                if not n.parent_node:
                    root_height = n.distance_from_root()
                    tree_height = root_height + n.distance_from_tip()
                    parent_height = 0
                else:
                    parent_height = n.parent_node.distance_from_root()
                node_height = round(node_height, 8)
                parent_height = round(parent_height, 8)

                if not node_height in blocks:
                    blocks[node_height] = []

                dists.append((n, parent_height, node_height))

            for time in blocks:
                for (node, parent_h, node_h) in dists:
                    if parent_h < time <= node_h:
                        blocks[time].append(node)

            dists.sort(key=lambda x: x[2])
            return (blocks, dists)

        def _weight_by_branches(blocks):
            intervals = sorted(blocks.keys())
            weighted_intervals = [0] + [None] * (len(intervals) - 1)
            for i in range(1, len(intervals)):
                time_range = intervals[i] - intervals[i - 1]
                num_branches = len(blocks[intervals[i]])
                weighted_range = time_range * num_branches
                weighted_intervals[i] = weighted_range + weighted_intervals[i
                        - 1]
            return weighted_intervals

        def _get_time(blocks, weights=None, verbose=False):
            d = sorted(blocks.keys())
            if weights:
                samp = random.uniform(weights[0], weights[-1])
                for i in range(len(weights) - 1):
                    if weights[i + 1] >= samp > weights[i]:
                        interval = weights[i + 1] - weights[i]
                        proportion = (samp - weights[i]) / interval
                        break
                drange = d[i + 1] - d[i]
                time = drange * proportion + d[i]
            else:
                time = random.uniform(d[0], d[-1])

            if verbose: 
                print 'LGT event at time: {0}'.format(time)

            return time

        def _choose_prune_and_regraft_nodes(time, blocks,
                disallow_sibling_SPRs, verbose=False):
            matching_branches = [x for x in dists if x[1] < time < x[2]]

            prune = random.sample(matching_branches, 1)[0]

            if disallow_sibling_SPRs:
                siblings = prune[0].sister_nodes()
                for br in matching_branches:
                    if br[0] in siblings:
                        matching_branches.remove(br)

            matching_branches.remove(prune)

            if matching_branches == []:
                if verbose: print 'No non-sibling branches available'
                return (None, None)

            regraft = random.sample(matching_branches, 1)[0]

            prune_taxa = [n.taxon.label for n in prune[0].leaf_iter()]
            regraft_taxa = [n.taxon.label for n in regraft[0].leaf_iter()]
            if verbose:
                print 'Donor group = {0}'.format(regraft_taxa)
                print 'Receiver group = {0}'.format(prune_taxa)
            return (prune, regraft)

        def _add_node(tree, time, regraft_node):
            parent_node = regraft_node[0].parent_node
            new_node = parent_node.add_child(dpy.Node(), edge_length=time
                    - regraft_node[1])
            tree.reindex_subcomponent_taxa()
            tree.update_splits()
            return new_node

        def _prunef(tree, node):
            tree.prune_subtree(node, update_splits=False,
                               delete_outdegree_one=True)

        def _regraftf(
            tree,
            time,
            target_node,
            child_node,
            ):

            target_node.add_child(child_node[0], edge_length=child_node[2]
                                  - time)
            return tree

        # MAIN
        tree = self.convert_to_dendropy_tree()
        tree.is_rooted = utils.dpy.check_rooted(self.tree.newick)

        (blocks, dists) = _get_blocks(tree)
        if not time:
            weights = _weight_by_branches(blocks)
            time = _get_time(blocks, weights, verbose=verbose)
        (p, r) = _choose_prune_and_regraft_nodes(time, dists,
                disallow_sibling_SPRs=disallow_sibling_SPRs, verbose=verbose)

        if (p, r) == (None, None):
            return self.spr(disallow_sibling_SPRs=disallow_sibling_SPRs)

        new_node = _add_node(tree, time, r)
        _prunef(tree, p[0])

        _prunef(tree, r[0])

        _regraftf(tree, time, new_node, p)

        _regraftf(tree, time, new_node, r)

        tree.reindex_subcomponent_taxa()
        tree.update_splits()

        newick = self.dendropy_as_newick(tree)
        if tree.is_rooted:
            newick = '[&R] ' + newick

        tree_copy = self.tree.copy() 
        tree_copy.newick = newick
        self.tree = tree_copy
        return self.tree


class Tree(object):

    """ Class for storing the results of phylogenetic inference """

    name_regex = re.compile('([A-Za-z0-9\-_]+).([A-Za-z0-9\-_]+)(?=_phyml_)')
    score_regex = re.compile('(?<=Log-likelihood: ).+')

    def __init__(
        self,
        newick=None,
        score=0,
        program=None,
        name=None,
        output=None,
        rooted=None,
        ):

        self.newick = newick
        self.score = score
        self.program = program
        self.name = name
        self.output = output
        self.rooted = utils.dpy.check_rooted(newick)

    def __repr__(self):
        return '{0}{1}'.format(self.__class__.__name__, (self.newick if self.newick else '(None)'))

    def __str__(self):
        """ Represents the object's information inside a newick comment, so is
        still interpretable by a (good) newick parser """

        s = '[Tree Object:\n'
        if self.name:
            s += 'Name:\t' + self.name + '\n'
        s += 'Program:\t{0}\n'.format(self.program) \
            + 'Score:\t{0}\n'.format(self.score) \
            + 'Rooted:\t{0}\n'.format(self.rooted) \
            + 'Tree:\t]{0}\n'.format(self.newick)
        return s

    def __eq__(self, other):
        equal = True
        if not self.name == other.name:
            return False
        if not self.newick == other.newick:
            return False
        if not self.program == other.program:
            return False
        if not self.score == other.score:
            return False
        if not self.rooted == other.rooted:
            return False
        if not self.output == other.output:
            return False
        return equal

    def copy(self):
        copy = self.__new__(type(self))
        copy.__dict__ = {key: value for (key, value) in self.__dict__.items()}
        return copy

    def read_from_file(self, infile, name=None):
        """ This and the write_to_file function allow the class to be easily
        stored and reconstituted without using a pickle or JSON """

        program = None
        tree = None
        score = None
        self.name = name
        reader = open(infile)
        try:
            for line in reader:
                line = [l.rstrip().replace(']', '') for l in line.split()]
                if not name and line[0] == 'Name:':
                    self.name = line[1]
                elif line[0] == 'Program:':
                    self.program = line[1]
                elif line[0] == 'Tree:':
                    self.newick = line[1]
                elif line[0] == 'Score:':
                    self.score = line[1]
        except IndexError:
            return
        return self

    def write_to_file(
        self,
        outfile,
        metadata=False,
        suppress_NHX=False,
        ):
        """ Writes a string representation of the object's contents to file.
        This can be read using read_from_file to reconstruct the Tree object, if
        metadata is included (i.e. metadata=True) """

        writer = open(outfile, 'w')
        if metadata:
            writer.write(str(self))
        else:

            writeable = self.newick
            if suppress_NHX:
                if writeable.startswith('[&R] '):
                    writeable = writeable[5:]
            if not writeable.endswith('\n'):
                writeable += '\n'
            writer.write(writeable)
        writer.close()
        return outfile

    def load_phyml_strings(
        self,
        tree,
        stats,
        name=None,
        program='phyml',
        ):

        score = float(self.score_regex.search(stats).group(1))
        self.program = program
        self.newick = tree
        self.output = stats
        self.score = score
        self.name = name
        self.rooted = utils.dpy.check_rooted(tree)

    def load_phyml_files(
        self,
        tree_file,
        stats_file,
        name=None,
        program='phyml',
        ):
        """ Loads phyml results into existing tree object - returns None """

        exit = False
        for f in (tree_file, stats_file):
            try:
                filecheck_and_raise(f)
            except FileError, e:
                print e
                exit = True

        if exit:
            print 'Results were not loaded'
            raise FileError()

        if not name:
            name = self.name_regex.search(tree_file).group(1)
        newick = open(tree_file).read()
        stats = open(stats_file).read()
        self.load_phyml_strings(newick, stats, name=name, program=program)

    @classmethod
    def new_tree_from_phyml_results(
        cls,
        tree_file,
        stats_file,
        program='phyml',
        ):
        """ classmethod version of load_phyml_files - returns a new Tree object
        """

        new_tree = cls()
        new_tree.load_phyml_files(tree_file, stats_file, program=program)
        return new_tree

    @classmethod
    def new_tree_from_phyml_strings(
        cls,
        tree,
        stats,
        program='phyml',
        ):

        new_tree = cls()
        new_tree.load_phyml_strings(tree, stats, program=program)
        return new_tree

    def scale(self, scale_factor):
        return TreeManip(self).scale(scale_factor)

    def strip(self):
        return TreeManip(self).strip()

    def print_plot(self):
        utils.dpy.print_plot(self)

    def rfdist(self, other):
        s = utils.dpy.convert_to_dendropy_tree(self)
        o = utils.dpy.convert_to_dendropy_tree(other)
        return utils.dpy.get_rf_distance(s, o)

    def wrfdist(self, other):
        s = utils.dpy.convert_to_dendropy_tree(self)
        o = utils.dpy.convert_to_dendropy_tree(other)
        return utils.dpy.get_wrf_distance(s, o)

    def eucdist(self, other):
        s = utils.dpy.convert_to_dendropy_tree(self)
        o = utils.dpy.convert_to_dendropy_tree(other)
        return utils.dpy.get_euc_distance(s, o)

    def geodist(self, other):
        gtp = GTP()
        return gtp.pairwise(self, other)

    @classmethod
    def new_yule(
        self,
        nspecies,
        names=None,
        cf=False,
        ):
        g = TreeGen(nspecies, names, cf=cf)
        return g.yule()

    @classmethod
    def new_coal(
        self,
        nspecies,
        names=None,
        cf=False,
        ):
        g = TreeGen(nspecies, names, cf=cf)
        return g.coal()

    @classmethod
    def new_rtree(
        self,
        nspecies,
        names=None,
        cf=False,
        ):
        g = TreeGen(nspecies, names, cf=cf)
        return g.rtree()

    def gene_tree(self, scale_to=None):
        g = TreeGen(template=self)
        return g.gene_tree(scale_to)['gene_tree']


class TreeGen(object):

    def __init__(
        self,
        nspecies=16,
        names=None,
        template=None,
        cf=False,
        ):

        self.nspecies = nspecies
        if cf:
            self.names = random.sample(cfnames, nspecies)
        else:
            self.names = names or ['Sp{0}'.format(i) for i in range(1, nspecies
                                   + 1)]
        if template and not isinstance(template, Tree):
            raise TypeError('template should be \'Tree\' object. Got',
                            type(tree))
        self.template = template

    def coal(self):
        taxon_set = dpy.TaxonSet(self.names)
        tree = dpy.treesim.pure_kingman(taxon_set)
        newick = '[&R] ' + utils.dpy.convert_dendropy_to_newick(tree)
        return Tree(newick)

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

        t = self.template or self.yule()
        tree = utils.dpy.convert_to_dendropy_tree(t)

        for leaf in tree.leaf_iter():
            leaf.num_genes = 1

        tree_height = tree.seed_node.distance_from_root() \
            + tree.seed_node.distance_from_tip()

        if scale_to:
            population_size = tree_height / scale_to

        for edge in tree.preorder_edge_iter():
            edge.pop_size = population_size

        gene_tree = dpy.treesim.constrained_kingman(tree)[0]

        if trim_names:
            for leaf in gene_tree.leaf_iter():
                leaf.taxon.label = leaf.taxon.label.replace('\'', '').split('_'
                        )[0]

        newick = '[&R] ' + utils.dpy.convert_dendropy_to_newick(gene_tree)

        return {'gene_tree': Tree(newick), 'species_tree': t}

    def rtree(self):
        m = TreeManip(self.yule())
        m.randomise_labels()
        return m.randomise_branch_lengths()

    def yule(self):
        taxon_set = dpy.TaxonSet(self.names)
        tree = dpy.treesim.uniform_pure_birth(taxon_set)
        newick = '[&R] ' + utils.dpy.convert_dendropy_to_newick(tree)
        return Tree(newick)


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
