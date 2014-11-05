#!/usr/bin/env python
from __future__ import print_function

# standard lib
import sys
import itertools
import random
import tempfile
import timeit

# third party
import numpy as np

# treeCl
from tree import Tree
from distance_matrix import DistanceMatrix
from alignment import Alignment
from utils import fileIO, flatten_list
from utils.decorators import lazyprop
from utils.printing import print_and_return
from errors import OptionError, optioncheck, directorycheck
from constants import SORT_KEY, PHYML_MEMORY_MIN


class NoRecordsError(Exception):
    def __init__(self, file_format, input_dir, compression):
        self.file_format = file_format
        self.input_dir = input_dir
        self.compression = compression

    def __str__(self):
        msg = ('No records were found in {0} matching\n'
               '\tfile_format = {1}\n'
               '\tcompression = {2}'.format(self.input_dir,
                                            self.file_format, self.compression))
        return msg


class Collection(object):
    """ Call:

    c = Collection(inut_dir, file_format, datatype, tmpdir ...)
    c.calc_distances(), c.calc_TC_trees(), ...
    dm = c.distance_matrix('geo')
    cl = Clustering(dm)
    k = cl.spectral(4, prune='estimate', local_scale=7)
    p = Partition(k) """

    def __init__(
            self,
            records=None,
            input_dir=None,
            trees_dir=None,
            file_format='fasta',
            compression=None,
    ):

        self._records = None
        self._input_files = None

        if records:
            self.records = records

        elif input_dir:
            directorycheck(input_dir)
            optioncheck(file_format, ['fasta', 'phylip'])
            self.records = self.read_alignments(input_dir,
                                                file_format,
                                                compression)

        else:
            raise Exception('Provide a list of records, '
                            'or the path to a set of alignments')

        if trees_dir:
            self.read_trees(trees_dir)

        if not self.records:
            raise NoRecordsError(file_format, input_dir, compression)


    def __len__(self):
        if hasattr(self, 'records'):
            return len(self.records)
        return 0

    def __getitem__(self, i):
        if hasattr(self, 'records'):
            return self.records[i]

    @property
    def records(self):
        """ Returns a list of records in SORT_KEY order """
        return [self._records[i] for i in range(len(self._records))]

    @records.setter
    def records(self, records):
        """ Sets a dictionary of records keyed by SORT_KEY order """
        self._records = dict(enumerate(records))

    @property
    def trees(self):
        """ Returns a list of trees in SORT_KEY order """
        try:
            return [rec.tree for rec in self]
        except ValueError:
            return []

    # @trees.setter
    # def trees(self, trees):
    # trees = dpy.TreeList(trees)
    # sorting_lambda = lambda x: SORT_KEY(x.name)
    # trees.sort(key=sorting_lambda)
    # for rec, tree in zip(self.records, trees):
    # assert rec.get_namespace() == tree.name
    #         rec.tree = tree

    def num_species(self):
        """ Returns the number of species found over all records
        """
        all_headers = reduce(lambda x, y: set(x) | set(y),
                             (rec.get_names() for rec in self.records))
        return len(all_headers)

    def read_alignments(self, input_dir, file_format, compression=None):
        """ Get list of alignment files from an input directory *.fa, *.fas and
        *.phy files only

        Stores in self.files """

        optioncheck(compression, [None, 'gz', 'bz2'])

        if file_format == 'fasta':
            extensions = ['fa', 'fas', 'fasta']

        elif file_format == 'phylip':
            extensions = ['phy']

        if compression:
            extensions = ['.'.join([x, compression]) for x in extensions]

        files = fileIO.glob_by_extensions(input_dir, extensions)
        files.sort(key=SORT_KEY)
        self._files = files
        records = []

        for f in files:
            if compression is not None:
                tmpdir = tempfile.mkdtemp()
                _, tmpfile = tempfile.mkstemp(dir=tmpdir)
                with fileIO.freader(f, compression) as reader:
                    with fileIO.fwriter(tmpfile) as writer:
                        for line in reader:
                            writer.write(line)
                try:
                    record = Alignment(tmpfile, file_format, True)
                except RuntimeError:
                    record = Alignment(tmpfile, file_format, False)
                finally:
                    os.remove(tmpfile)
                    os.rmdir(tmpdir)
            else:
                try:
                    record = Alignment(f, file_format, True)
                except RuntimeError:
                    record = Alignment(f, file_format, False)

            record.name = (fileIO.strip_extensions(f))
            record.fast_compute_distances()
            records.append(record)

        return records

    def read_trees(self, input_dir):
        """ Read a directory full of tree files, matching them up to the
        already loaded alignments """

        extensions = ['nwk', 'tree']
        files = fileIO.glob_by_extensions(input_dir, extensions)
        trees = [Tree.read_from_file(file_) for file_ in files]

    def fast_calc_distances(self):
        """ Calculates within-alignment pairwise distances and variances for every
        alignment. Uses fast Jukes-Cantor method.
        :return: void"""
        for rec in self.records:
            print_and_return("Calculating fast distances for {}".format(rec.name))
            rec.fast_compute_distances()

    def calc_distances(self):
        """ Calculates within-alignment pairwise distances and variances for every
        alignment. Uses slow ML optimisation, and depends on the Alignment record
        having had appropriate ML parametric models set up in advance.
        :return: void
        """
        for rec in self.records:
            print_and_return("Calculating distances for {}".format(rec.name))
            rec.compute_distances()

    def calc_pll_trees(self, threads=1):
        """ Use pllpy to calculate maximum-likelihood trees
        :return: void
        """
        for rec in self.records:
            print_and_return("Calculating ML tree for {}".format(rec.name))
            if rec.is_dna():
                model = 'GTR'
            else:
                model = 'LGX'
            result = rec.pll_optimise('{}, {} = 1 - {}'.format(model, rec.name, len(rec)), rec.tree.newick, nthreads=threads)
            freqs = result['partitions'][0]['frequencies']
            tree = result['tree']
            alpha = result['partitions'][0]['alpha']
            rec.set_substitution_model('GTR' if rec.is_dna() else 'LG08')
            rec.set_gamma_rate_model(4, alpha)
            rec.set_frequencies(freqs)
            if rec.is_dna():
                record.set_rates(result['partitions'][0]['rates'], 'ACGT')
            rec.initialise_likelihood(tree)

    def get_tree_distance_matrix(self, metric, **kwargs):
        """ Generate a distance matrix from a fully-populated Collection """
        return DistanceMatrix(self.trees, metric, **kwargs)

    def permuted_copy(self):
        """ Return a copy of the collection with all alignment columns permuted
        """

        def take(n, iterable):
            return [iterable.next() for _ in range(n)]

        def items_subset(keys, d):
            return [(k, d[k]) for k in keys]

        concat = Concatenation(self, range(len(self)))
        sites = concat.alignment.get_sites()
        random.shuffle(sites)
        d = dict(zip(concat.alignment.get_names(), [iter(x) for x in zip(*sites)]))

        new_seqs = []
        for l in concat.lengths:
            new_seqs.append(dict([(k, ''.join(take(l, d[k]))) for k in d]))

        records = []
        for (k, d) in zip(concat.headers, new_seqs):
            records.append(items_subset(k, d))

        permutation = self.__class__(
            records=[Alignment(seqs, dtype) for (seqs, dtype) in zip(records, concat.datatypes)])
        for rec, name in zip(permutation, concat.names):
            rec.name = name

        return permutation


class Concatenation(object):
    """docstring for Concatenation"""

    def __init__(self, collection, indices):
        super(Concatenation, self).__init__()
        if any((x > len(collection)) for x in indices):
            raise ValueError('Index out of bounds in {}'.format(indices))
        if any((x < 0) for x in indices) < 0:
            raise ValueError('Index out of bounds in {}'.format(indices))
        if any((not isinstance(x, int)) for x in indices):
            raise ValueError('Integers only in indices, please: {}'
                             .format(indices))
        self.collection = collection
        self.indices = sorted(indices)

    @lazyprop
    def distances(self):
        return [self.collection.records[i].get_distances() for i in self.indices]

    @lazyprop
    def datatypes(self):
        return ['dna' if self.collection.records[i].is_dna() else 'protein' for i in self.indices]

    @lazyprop
    def alignment(self):
        return Alignment([self.collection[i] for i in self.indices])

    @lazyprop
    def names(self):
        return [self.collection.records[i].name for i in self.indices]

    @lazyprop
    def lengths(self):
        return [len(self.collection.records[i]) for i in self.indices]

    @lazyprop
    def headers(self):
        return [self.collection.records[i].get_names() for i in self.indices]

    @lazyprop
    def coverage(self):
        total = float(self.collection.num_species())
        return [len(self.collection.records[i]) / total for i in self.indices]

    @lazyprop
    def trees(self):
        return [self.collection.records[i].tree for i in self.indices]

    @lazyprop
    def mrp_tree(self):
        trees = [tree.newick for tree in self.trees]
        return Tree(Alignment().get_mrp_supertree(trees))

    def _get_tree_collection_strings(self, scale=1):
        """ Function to get input strings for tree_collection
        tree_collection needs distvar, genome_map and labels -
        these are returned in the order above
        """

        # aliases
        num_matrices = len(self.distances)
        label_set = reduce(lambda x, y: x.union(y), (set(l) for l in self.headers))
        labels_len = len(label_set)

        # labels string can be built straight away
        labels_string = '{0}\n{1}\n'.format(labels_len, ' '.join(label_set))

        # distvar and genome_map need to be built up
        distvar_list = [str(num_matrices)]
        genome_map_list = ['{0} {1}'.format(num_matrices, labels_len)]

        # build up lists to turn into strings
        for i in range(num_matrices):
            labels = self.headers[i]
            dim = len(labels)
            matrix = self.distances[i].copy()
            if scale:
                matrix[np.triu_indices(dim, 1)] *= scale
                matrix[np.tril_indices(dim, -1)] *= scale * scale

            if isinstance(matrix, np.ndarray):
                matrix_string = '\n'.join([' '.join(str(x) for x in row)
                                           for row in matrix]) + '\n'
            else:
                matrix_string = matrix
            distvar_list.append('{0} {0} {1}\n{2}'.format(dim, i + 1,
                                                          matrix_string))
            genome_map_entry = ' '.join((str(labels.index(lab) + 1)
                                         if lab in labels else '-1')
                                        for lab in label_set)
            genome_map_list.append(genome_map_entry)

        distvar_string = '\n'.join(distvar_list)
        genome_map_string = '\n'.join(genome_map_list)

        return distvar_string, genome_map_string, labels_string

    def tree_collection(self,
                        niters=5,
                        keep_topology=False,
                        quiet=True,
                        guide_tree=None,
                        scale=1):

        import tree_collection

        guide_tree = guide_tree.copy()
        for e in guide_tree.postorder_edge_iter():
            if e.length is None:
                if e.head_node == guide_tree.seed_node:
                    e.length = 0.0
                else:
                    e.length = 1.0

        if not guide_tree.is_rooted:
            guide_tree.reroot_at_midpoint()
        if not guide_tree.is_rooted:
            raise Exception('Couldn\'t root the guide tree')

        dv, gm, lab = self._get_tree_collection_strings(scale)
        output_tree, score = tree_collection.compute(dv, gm, lab, guide_tree.scale(scale).newick,
                                                     niters, keep_topology,
                                                     quiet)

        return Tree(output_tree), score


    def qfile(self, dna_model='DNA', protein_model='LG', sep_codon_pos=False,
              ml_freqs=False, eq_freqs=False):
        from_ = 1
        to_ = 0
        qs = list()
        if ml_freqs:
            dna_model += 'X'
            protein_model += 'X'
        if eq_freqs and not ml_freqs:
            protein_model += 'F'

        models = dict(dna=dna_model, protein=protein_model)
        for length, name, datatype in zip(self.lengths, self.names,
                                          self.datatypes):
            to_ += length
            if datatype == 'dna' and sep_codon_pos:
                qs.append('{}, {} = {}-{}/3'.format(models[datatype], name, from_,
                                                    to_))
                qs.append('{}, {} = {}-{}/3'.format(models[datatype], name, from_ + 1,
                                                    to_))
                qs.append('{}, {} = {}-{}/3'.format(models[datatype], name, from_ + 2,
                                                    to_))
            else:
                qs.append('{}, {} = {}-{}'.format(models[datatype], name, from_,
                                                  to_))
            from_ += length
        return '\n'.join(qs)

    def pll_optimise(self, partitions, tree=None, model=None, nthreads=1, **kwargs):
        if tree is None:
            tree = self.mrp_tree.newick
        return self.alignment.pll_optimise(partitions, tree, model, nthreads, **kwargs)

    def paml_partitions(self):
        return 'G {} {}'.format(len(self.lengths),
                                ' '.join(str(x) for x in self.lengths))


class Scorer(object):
    """ Takes an index list, generates a concatenated SequenceRecord, calculates
    a tree and score """

    def __init__(
            self,
            collection,
            verbosity=0,
            populate_cache=True,
    ):

        self.collection = collection
        self.verbosity = verbosity
        self.cache = {}
        self.history = []
        if populate_cache:
            self.populate_cache()

    @property
    def records(self):
        return self.collection.records

    def add_partition_list(self, partition_list):
        """ Calculates concatenated trees for a list of Partitions """
        index_tuples = list(itertools.chain(*[partition.get_membership()
                                              for partition in partition_list]))
        missing = sorted(set(index_tuples).difference(self.cache.keys()))
        self._add_index_tuple_list(missing)

    def _add_index_tuple_list(self, index_tuple_list):
        if self.lsf and not self.analysis == 'TreeCollection':
            supermatrices = [self.concatenate(index_tuple).alignment
                             for index_tuple in index_tuple_list]

            trees = runLSFPhyml(supermatrices,
                                self.tmpdir,
                                analysis=self.analysis,
                                verbosity=self.verbosity,
                                strategy='dynamic',
                                minmem=PHYML_MEMORY_MIN,
                                debug=self.debug,
            )
            for index_tuple, tree in zip(index_tuple_list, trees):
                self.cache[index_tuple] = tree
        else:
            for index_tuple in index_tuple_list:
                self.add(index_tuple)

    def add(self, index_tuple):
        """ Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree """

        if index_tuple in self.cache:
            return self.cache[index_tuple]

        if len(index_tuple) == 1:
            alignment = self.collection[index_tuple[0]]
        else:
            concat = self.concatenate(index_tuple)
            alignment = concat.alignment

        if self.analysis == 'TreeCollection':
            guidetrees = [self.records[n].tree for n in
                          index_tuple][:self.max_guidetrees]
            tree = alignment.tree_collection(
                guide_tree=guidetrees[0])

        else:
            tree = runPhyml(alignment,
                            self.tmpdir,
                            analysis=self.analysis,
                            verbosity=self.verbosity, )
            if self.verbosity == 1:
                print()

        self.cache[index_tuple] = tree
        return tree

    def concatenate(self, index_tuple):
        """ Returns a Concatenation object that stitches together
        the alignments picked out by the index tuple """
        return Concatenation(self.collection, index_tuple)

    def update_history(self, score, index_tuple):
        """ Used for logging the optimiser """
        time = timeit.default_timer()
        self.history.append([time, score, index_tuple, len(index_tuple)])

    def print_history(self, fh=sys.stdout):
        """ Used for logging the optimiser """
        for iteration, (time, score, index_tuple, nclusters) in enumerate(
                self.history):
            fh.write(str(iteration) + "\t")
            fh.write(str(time) + "\t")
            fh.write(str(score) + "\t")
            fh.write(str(index_tuple) + "\t")
            fh.write(str(nclusters) + "\n")

    def clear_history(self):
        """ Used for logging the optimiser: clears the log """
        self.history = []

    def members(self, index_tuple):
        """ Gets records by their index, contained in the index_tuple """
        return [self.records[n] for n in index_tuple]

    def populate_cache(self):
        """ Adds all single-record trees to the cache """
        to_calc = []
        for i, rec in enumerate(self.records):
            key = (i,)
            if rec.tree is None:
                to_calc.append(key)
                continue
            if rec.tree.program.startswith('phyml+'):
                analysis = rec.tree.program[6:]
            else:
                analysis = rec.tree.program
            if analysis == self.analysis:
                tree = rec.tree
                self.cache[key] = tree
            else:
                to_calc.append(key)
        self._add_index_tuple_list(to_calc)


    def score(self, partition_object, history=True, **kwargs):
        """ Generates the index lists of the Partition object, gets the score
        for each one, and returns the sum """

        inds = partition_object.get_membership()
        self.add_partition_list([partition_object])
        likelihood = sum([self.add(index_tuple, **kwargs).score
                          for index_tuple in inds])
        if history is True:
            self.update_history(likelihood, inds)
        return likelihood

    def simulate(self, index_tuple, model=None, lsf=False, ntimes=1):
        """ Simulate a group of sequence alignments using ALF. Uses one of
        {(GCB, JTT, LG, WAG - protein), (CPAM, ECM and ECMu - DNA)}, WAG by
        default. TO DO: add parameterised models when I have a robust (probably
        PAML) method of estimating them from alignment+tree """

        if self.datatype == 'protein':  # set some defaults
            model = model or 'WAG'
            optioncheck(model, [
                'CPAM',
                'ECM',
                'ECMu',
                'WAG',
                'JTT',
                'GCB',
                'LG',
            ])
        else:
            model = model or 'GTR'
            try:
                optioncheck(model, ['CPAM', 'ECM', 'ECMu', 'GTR'])
            except OptionError, e:
                print('Choose a DNA-friendly model for simulation:\n', e)
                return

        member_records = self.members(index_tuple)
        concat = self.concatenate(index_tuple).alignment
        (lengths, names) = zip(*[(rec.seqlength, rec.name) for rec in
                                 member_records])
        full_length = sum(lengths)
        concat.tree = self.add(index_tuple)

        if lsf and ntimes > 1:
            simulated_records = lsf_simulate_from_record(
                concat,
                ntimes,
                length=full_length,
                tmpdir=self.tmpdir,
                model=model,
                split_lengths=lengths,
                gene_names=names,
            )

        else:
            simulated_records = simulate_from_record(
                concat,
                length=full_length,
                tmpdir=self.tmpdir,
                model=model,
                split_lengths=lengths,
                gene_names=names,
            )

        return simulated_records

    def simulate_from_result(self,
                             partition_object, lsf=False,
                             ntimes=1, **kwargs
    ):
        """ Simulates a set of records using parameters estimated when
        calculating concatenated trees from the Partition object """
        inds = partition_object.get_membership()

        if lsf and ntimes > 1:
            multiple_results = [self.simulate(ind, lsf=lsf, ntimes=ntimes)
                                for ind in inds]
            return [flatten_list(result)
                    for result in zip(*multiple_results)]

        else:
            return [flatten_list([self.simulate(ind, **kwargs)
                                  for ind in inds])
                    for _ in range(ntimes)]


    def dump_cache(self, filename):
        with open(filename, 'w') as outfile:
            for k, v in self.cache.items():
                outfile.write('{}\t{}\n'.format(k, str(v)))

    def load_cache(self, filename):
        d = {}
        with open(filename, 'r') as infile:
            for line in infile:
                k, s = line.rstrip().split('\t')
                t = Tree.gen_from_text(s)
                d[k] = t
        self.cache = d
