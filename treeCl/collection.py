#!/usr/bin/env python
from __future__ import print_function

# standard lib
import os
import sys
import random
import tempfile
import timeit

# third party
import numpy as np

# treeCl
from tree import Tree
from distance_matrix import DistanceMatrix
from alignment import Alignment
from utils import fileIO, setup_progressbar
from utils.decorators import lazyprop
from utils.printing import print_and_return
from errors import optioncheck, directorycheck
from constants import SORT_KEY, PLL_RANDOM_SEED


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

    c = Collection(input_dir, file_format, datatype, tmpdir ...)
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
    # rec.tree = tree

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

        else:
            extensions = []

        if compression:
            extensions = ['.'.join([x, compression]) for x in extensions]

        files = fileIO.glob_by_extensions(input_dir, extensions)
        files.sort(key=SORT_KEY)
        self._input_files = files
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
        pbar = setup_progressbar('Calculating distances', len(self))
        pbar.start()
        for i, rec in enumerate(self.records):
            rec.compute_distances()
            pbar.update(i)
        pbar.finish()

    def calc_pll_trees(self, threads=1):
        """ Use pllpy to calculate maximum-likelihood trees
        :return: void
        """
        pbar = setup_progressbar('Calculating ML trees', len(self))
        pbar.start()
        for i, rec in enumerate(self.records):
            if rec.is_dna():
                model = 'GTR'
            else:
                model = 'LGX'
            result = rec.pll_optimise('{}, {} = 1 - {}'.format(model, rec.name, len(rec)), rec.tree.newick,
                                                               nthreads=threads, seed=PLL_RANDOM_SEED)
            freqs = result['partitions'][0]['frequencies']
            tree = result['tree']
            alpha = result['partitions'][0]['alpha']
            rec.set_substitution_model('GTR' if rec.is_dna() else 'LG08')
            rec.set_gamma_rate_model(4, alpha)
            rec.set_frequencies(freqs)
            if rec.is_dna():
                rec.set_rates(result['partitions'][0]['rates'], 'ACGT')
            rec.initialise_likelihood(tree)
            pbar.update(i)
        pbar.finish()

    def get_tree_distance_matrix(self, metric, **kwargs):
        """ Generate a distance matrix from a fully-populated Collection """
        return DistanceMatrix(self.trees, metric, **kwargs)

    def permuted_copy(self):
        """ Return a copy of the collection with all alignment columns permuted
        """

        def take(n, iterable):
            return [iterable.next() for _ in range(n)]

        def items_subset(kys, dct):
            return [(ky, dct[ky]) for ky in kys]

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
        al = Alignment([self.collection[i] for i in self.indices])
        al.fast_compute_distances()
        return al

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

        guide_tree = self.alignment.tree

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
        tree_string = guide_tree.scale(scale).newick

        return distvar_string, genome_map_string, labels_string, tree_string

    def minsq_tree(self,
                   niters=5,
                   keep_topology=False,
                   quiet=True,
                   scale=1):

        dv, gm, lab, tree_string = self._get_tree_collection_strings(scale)

        import tree_collection
        output_tree, score = tree_collection.compute(dv, gm, lab, tree_string,
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
            tree = self.alignment.tree.newick
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
    ):

        self.collection = collection
        self.verbosity = verbosity
        self.minsq_cache = {}
        self.lnl_cache = {}
        self.history = []

    @property
    def records(self):
        return self.collection.records

    def get_minsq_partition(self, partition_list):
        """ Calculates concatenated trees for a list of Partitions """
        index_tuples = partition.get_membership()
        return self._get_minsq_index_tuple_list(index_tuples)


    def get_lnl_partition(self, partition):
        """ Calculates concatenated trees for a list of Partitions """
        index_tuples = partition.get_membership()
        return self._get_lnl_index_tuple_list(index_tuples)

    def _get_lnl_index_tuple_list(self, index_tuple_list):
        """
        Does maximum-likelihood tree optimisation for the alignments specified
        by the tuple list.
        :param index_tuple_list:
        :return:
        """
        return [self._get_lnl(index_tuple) for index_tuple in index_tuple_list]

    def _get_minsq_index_tuple_list(self, index_tuple_list):
        """
        Does tree estimation by minimum squared distance for the alignments
        specified by the tuple list.
        :param index_tuple_list:
        :return:
        """
        return [self._get_minsq(index_tuple) for index_tuple in index_tuple_list]

    def _get_lnl(self, index_tuple):
        """
        Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree
        :param index_tuple: tuple of indexes at which to find the alignments.
        :return:
        """
        try:
            return self.lnl_cache[index_tuple]
        except KeyError:
            conc = self.concatenate(index_tuple)
            partitions = conc.qfile(dna_model="GTR", protein_model="LGX")
            result = self.lnl_cache[index_tuple] = conc.pll_optimise(partitions)
            return result

    def _get_minsq(self, index_tuple):
        """
        Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree
        :param index_tuple: tuple of indexes at which to find the alignments.
        :return:
        """
        try:
            return self.minsq_cache[index_tuple]
        except KeyError:
            conc = self.concatenate(index_tuple)
            result = self.minsq_cache[index_tuple] = conc.minsq_tree()
            return result

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

    def score(self, partition, criterion):
        """
        Return the score for a partition - either the sum of log likelihoods,
        or the total min squares dimensionless fit index
        :param partition: Partition object
        :param criterion: either 'minsq' or 'lnl'
        :return: score (float)
        """
        optioncheck(criterion, ['lnl', 'minsq'])
        results = (self.get_lnl_partition(partition) if criterion == 'lnl'
                   else self.get_minsq_partition(partition))
        return results  # TODO: how to sum scores from this?
