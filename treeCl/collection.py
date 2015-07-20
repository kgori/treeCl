#!/usr/bin/env python
from __future__ import print_function

# standard lib
import glob
import itertools
import json
import math
import os
import random
import sys

# third party
import numpy as np
from scipy.spatial.distance import squareform

# treeCl
from .alignment import Alignment
from .concatenation import Concatenation
from .constants import SORT_KEY, PLL_RANDOM_SEED
from .distance_matrix import DistanceMatrix
from .errors import optioncheck, directorycheck
from . import tasks
from .parameters import PartitionParameters
from .partition import Partition
from .parutils import get_client, parallel_map, sequential_map
from .tree import Tree
from .utils import fileIO, setup_progressbar, model_translate
from .utils.decorators import lazyprop

import json
import sys
import logging
logger = logging.getLogger(__name__)

def gapmask(simseqs, origseqs):
    """
    :param sims: list of (header, sequence) tuples of simulated sequences [no gaps]
    :param aln: list of (header, sequence) tuples of original sequences
    :return:
    """
    import numpy as np
    simdict = dict(simseqs)
    origdict = dict(origseqs)
    for k in origdict:
        origseq = np.array(list(origdict[k]))
        gap_pos = np.where(origseq=='-')
        simseq = np.array(list(simdict[k]))
        simseq[gap_pos] = '-'
        simdict[k] = ''.join(simseq)
    return list(simdict.items())


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
            param_dir=None,
            file_format='fasta',
            compression=None,
            header_grep=None,
            show_progressbars=True,
    ):

        self._records = None
        self._input_files = None
        self.show_progressbars=show_progressbars

        if records is not None:
            self.records = records

        elif input_dir is not None:
            input_dir = os.path.abspath(input_dir)
            directorycheck(input_dir)
            optioncheck(file_format, ['fasta', 'phylip'])
            self.records = self.read_alignments(input_dir,
                                                file_format,
                                                header_grep,
                                                compression)

        else:
            raise Exception('Provide a list of records, '
                            'or the path to a set of alignments')

        if param_dir is not None:
            self.read_parameters(param_dir)

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

    @lazyprop
    def trees(self):
        """ Returns a list of trees in SORT_KEY order """
        try:
            return [rec.tree for rec in self]
        except ValueError:
            return []

    @lazyprop
    def names(self):
        """
        Returns a list of sequence record names in SORT_KEY order
        """
        try:
            return [rec.name for rec in self]
        except ValueError:
            return []

    @lazyprop
    def distances(self):
        try:
            return [rec.parameters.partitions.distances for rec in self]
        except ValueError:
            return []

    @lazyprop
    def variances(self):
        try:
            return [rec.parameters.partitions.variances for rec in self]
        except ValueError:
            return []

    @lazyprop
    def frequencies(self):
        try:
            return [rec.parameters.partitions.frequencies for rec in self]
        except ValueError:
            return []

    @lazyprop
    def alphas(self):
        try:
            return [rec.parameters.partitions.alpha for rec in self]
        except ValueError:
            return []

    @lazyprop
    def datatypes(self):
        try:
            return ['dna' if rec.is_dna() else 'protein' for rec in self]
        except ValueError:
            return []

    @lazyprop
    def lengths(self):
        try:
            return [len(rec) for rec in self]
        except ValueError:
            return []

    @lazyprop
    def headers(self):
        try:
            return [rec.get_names() for rec in self]
        except ValueError:
            return []

    @lazyprop
    def mrp_tree(self):
        trees = [tree.newick if hasattr('newick', tree) else tree for tree in self.trees]
        return Alignment().get_mrp_supertree(trees)

    def num_species(self):
        """ Returns the number of species found over all records
        """
        all_headers = reduce(lambda x, y: set(x) | set(y),
                             (rec.get_names() for rec in self.records))
        return len(all_headers)

    def read_alignments(self, input_dir, file_format, header_grep=None, compression=None):
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

        if self.show_progressbars:
            pbar = setup_progressbar("Loading files", len(files), simple_progress=True)
            pbar.start()

        for i, f in enumerate(files):
            if compression is not None:
                with fileIO.TempFile() as tmpfile:
                    with fileIO.freader(f, compression) as reader, fileIO.fwriter(tmpfile) as writer:
                        for line in reader:
                            writer.write(line)
                    try:
                        record = Alignment(tmpfile, file_format, True)
                    except RuntimeError:
                        record = Alignment(tmpfile, file_format, False)

            else:
                try:
                    record = Alignment(f, file_format, True)
                except RuntimeError:
                    record = Alignment(f, file_format, False)

            if header_grep:
                try:
                    datatype = 'dna' if record.is_dna() else 'protein'

                    record = Alignment([(header_grep(x), y) for (x, y) in record.get_sequences()], datatype)

                except TypeError:
                    raise TypeError("Couldn't apply header_grep to header\n"
                                    "alignment number={}, name={}\n"
                                    "header_grep={}".format(i, fileIO.strip_extensions(f), header_grep))
                except RuntimeError:
                    print('RuntimeError occurred processing alignment number={}, name={}'
                          .format(i, fileIO.strip_extensions(f)))
                    raise

            record.name = (fileIO.strip_extensions(f))
            records.append(record)
            if self.show_progressbars:
                pbar.update(i)
        if self.show_progressbars:
            pbar.finish()
        return records

    def read_parameters(self, input_dir):
        """ Read a directory full of tree files, matching them up to the
        already loaded alignments """

        if self.show_progressbars:
            pbar = setup_progressbar("Loading parameters", len(self.records))
            pbar.start()
        for i, rec in enumerate(self.records):
            hook = os.path.join(input_dir, '{}.json*'.format(rec.name))
            filename = glob.glob(hook)
            try:
                with fileIO.freader(filename[0]) as infile:
                    d = json.load(infile, parse_int=True)

                rec.parameters.construct_from_dict(d)

            except IOError, IndexError:
                continue

            finally:
                if self.show_progressbars:
                    pbar.update(i)
        if self.show_progressbars:
            pbar.finish()

    def write_parameters(self, output_dir, gz=False):
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except IOError as err:
                sys.stderr.write(err.message)
                raise err

        for rec in self.records:
            with fileIO.fwriter(os.path.join(output_dir, '{}.json'.format(rec.name)), gz=True) as outfile:
                rec.parameters.write(outfile, indent=4)

    def permuted_copy(self, partition=None):
        """ Return a copy of the collection with all alignment columns permuted
        """
        def take(n, iterable):
            return [iterable.next() for _ in range(n)]

        if partition is None:
            partition = Partition([1] * len(self))

        index_tuples = partition.get_membership()

        alignments = []
        for ix in index_tuples:
            concat = Concatenation(self, ix)
            sites = concat.alignment.get_sites()
            random.shuffle(sites)
            d = dict(zip(concat.alignment.get_names(), [iter(x) for x in zip(*sites)]))
            new_seqs = [[(k, ''.join(take(l, d[k]))) for k in d] for l in concat.lengths]

            for seqs, datatype, name in zip(new_seqs, concat.datatypes, concat.names):
                alignment = Alignment(seqs, datatype)
                alignment.name = name
                alignments.append(alignment)

        return self.__class__(records=sorted(alignments, key=lambda x: SORT_KEY(x.name)))

    def fast_calc_distances(self, batchsize=1, background=False):
        """
        Calculate fast approximate intra-alignment pairwise distances and variances using
        Jukes-Cantor closed formulae.
        :return: None (all side effects)
        """
        # Assemble argument lists
        args = []
        to_delete = []
        for rec in self:
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            args.append((filename,))

        # Dispatch work (either sequentially or in parallel)
        msg = 'Calculating fast distances'
        with fileIO.TempFileList(to_delete):
            client = get_client()
            if client is None:
                map_result = sequential_map(tasks.fast_calc_distances_task, args, msg)
            else:
                map_result = parallel_map(client, tasks.fast_calc_distances_task, args, msg, batchsize, background)
                if background:
                    return map_result

        # Process results
        pbar = setup_progressbar('Processing results', len(map_result))
        j = 0
        pbar.start()
        for i, result in enumerate(map_result):
            rec = self[i]
            distances = result['distances']
            variances = result['variances']
            tree = result['tree']
            rec.parameters.nj_tree = tree
            params = rec.parameters.partitions
            if params is None:
                params = PartitionParameters()
                rec.parameters.partitions = [params]
            params.distances = distances
            params.variances = variances
            pbar.update(i)
        pbar.finish()

    def calc_distances(self, batchsize=1, background=False):
        """
        Calculate fast approximate intra-alignment pairwise distances and variances using
        ML (requires ML models to have been set up using `calc_trees`).
        :return: None (all side effects)
        """
        # Assemble argument lists
        args = []
        to_delete = []
        for rec in self:
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            # Get input dict
            model = {'partitions': {}}
            data = {'alpha': rec.parameters.partitions.alpha, 'frequencies': rec.parameters.partitions.frequencies}
            if rec.is_dna():
                data['rates'] = rec.parameters.partitions.rates
            model['partitions'][0] = data
            args.append((model, filename))

        # Dispatch
        msg = 'Calculating ML distances'
        client = get_client()
        if client is None:
            map_result = sequential_map(tasks.calc_distances_task, args, msg)
        else:
            map_result = parallel_map(client, tasks.calc_distances_task, args, msg, batchsize, background)
            if background:
                return map_result

        # Process results
        with fileIO.TempFileList(to_delete):
            pbar = setup_progressbar('Processing results', len(map_result))
            j = 0
            pbar.start()
            for i, result in enumerate(map_result):
                rec = self[i]
                rec.parameters.partitions.distances = result['partitions'][0]['distances']
                rec.parameters.partitions.variances = result['partitions'][0]['variances']
                rec.parameters.nj_tree = result['nj_tree']
                pbar.update(j+1)
                j += 1
            pbar.finish()

    def fast_calc_trees(self, indices=None, batchsize=1, background=False):
        """
        Use FastTree to calculate maximum-likelihood trees
        :return: None (all side effects)
        """
        # Assemble argument lists
        if indices is None:
            indices = list(range(len(self)))
        args = []
        to_delete = []
        for i in indices:
            rec = self[i]
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            
            curr_args = (filename, rec.is_dna())
            args.append(curr_args)

        # Dispatch work
        msg = 'Calculating FastTree trees'
        client = get_client()
        if client is None:
            map_result = sequential_map(tasks.fasttree_task, args, msg)
        else:
            map_result = parallel_map(client, tasks.fasttree_task, args, msg, batchsize, background)
            if background:
                return map_result

        # Process results
        with fileIO.TempFileList(to_delete):
            pbar = setup_progressbar('Processing results', len(map_result))
            j = 0
            pbar.start()
            for i, result in zip(indices, map_result):
                rec = self[i]
                rec.parameters.construct_from_dict(result)
                pbar.update(j+1)
                j += 1
            pbar.finish()

    def calc_trees(self, model=None, threads=1, indices=None, tree_search=True, batchsize=1, output_dir=None, background=False):
        """
        Use pllpy to calculate maximum-likelihood trees
        :return: None (all side effects)
        """
        # Assemble argument lists
        if indices is None:
            indices = list(range(len(self)))
        args = []
        to_delete = []
        for i in indices:
            rec = self[i]
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            if model is None:
                model = ('DNA' if rec.is_dna() else 'LGX')
            if model == 'AUTOX':
                model = 'AUTO'
            partition = '{}, {} = 1 - {}'.format(model, rec.name, len(rec))
            tree = rec.parameters.nj_tree if rec.parameters.nj_tree is not None else True
            if output_dir is not None and os.path.isdir(output_dir):
                output_file = os.path.join(output_dir, '{}.json'.format(rec.name))
                curr_args = (filename, partition, tree, tree_search, threads, PLL_RANDOM_SEED, None, output_file)
            else:
                curr_args = (filename, partition, tree, tree_search, threads, PLL_RANDOM_SEED)
            args.append(curr_args)

        # Dispatch work
        msg = 'Calculating ML trees'
        client = get_client()
        if client is None:
            map_result = sequential_map(tasks.pll_task, args, msg)
        else:
            map_result = parallel_map(client, tasks.pll_task, args, msg, batchsize, background)
            if background:
                return map_result

        # Process results
        with fileIO.TempFileList(to_delete):
            pbar = setup_progressbar('Processing results', len(map_result))
            j = 0
            pbar.start()
            for i, result in zip(indices, map_result):
                rec = self[i]
                rec.parameters.construct_from_dict(result)
                pbar.update(j+1)
                j += 1
            pbar.finish()

    def concatenate(self, indices):
        return Concatenation(self, indices)

    def get_inter_tree_distances(self, metric, normalise=False, batchsize=100, background=False):
        """ Generate a distance matrix from a fully-populated Collection """
        array = _get_inter_tree_distances(metric, self.trees, normalise, batchsize, background)
        if background:  # return IPython.parallel map result object to the user before jobs are finished
            return array
        return DistanceMatrix.from_array(array, self.names)

    def get_tree_collection_strings(self, indices, scale=1, guide_tree=None):
        """ Function to get input strings for tree_collection
        tree_collection needs distvar, genome_map and labels -
        these are returned in the order above
        """

        # local lists
        distances = []
        variances = []
        headers = []
        for i in indices:
            distances.append(self.distances[i])
            variances.append(self.variances[i])
            headers.append(self.headers[i])

        num_matrices = len(indices)
        label_set = reduce(lambda x, y: x.union(y), (set(l) for l in headers))
        labels_len = len(label_set)

        # labels string can be built straight away
        labels_string = '{0}\n{1}\n'.format(labels_len, ' '.join(label_set))

        # distvar and genome_map need to be built up
        distvar_list = [str(num_matrices)]
        genome_map_list = ['{0} {1}'.format(num_matrices, labels_len)]

        # build up lists to turn into strings
        for i in range(num_matrices):
            labels = headers[i]
            dim = len(labels)
            dmatrix = np.array(distances[i])
            vmatrix = np.array(variances[i])
            matrix = np.zeros(dmatrix.shape)
            matrix[np.triu_indices(len(dmatrix), 1)] = dmatrix[np.triu_indices(len(dmatrix), 1)]
            matrix[np.tril_indices(len(vmatrix), -1)] = vmatrix[np.tril_indices(len(vmatrix), -1)]
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

        if guide_tree is None:
            guide_tree = Tree.new_iterative_rtree(labels_len, names=label_set, rooted=True)

        tree_string = guide_tree.scale(scale).newick

        return distvar_string, genome_map_string, labels_string, tree_string


def _get_inter_tree_distances(metric, trees, normalise=False, batchsize=100, background=False):
    # Assemble argument lists
    args = [(t1, t2, normalise) for (t1, t2) in itertools.combinations(trees, 2)]

    # Get task
    tasks_dict = dict(zip(['euc', 'geo', 'rf', 'wrf'],
                          [tasks.eucdist_task, tasks.geodist_task, tasks.rfdist_task, tasks.wrfdist_task]))
    task = tasks_dict[metric]

    # Dispatch
    msg = 'Inter-tree distances ({})'.format(metric)
    client = get_client()
    if client is None:
        map_result = sequential_map(task, args, msg)
    else:
        map_result = parallel_map(client, task, args, msg, batchsize, background)
        if background:
            return map_result
        map_result = list(map_result)

    return squareform(map_result)


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

    @property
    def records(self):
        return self.collection.records

    def add_lnl_partitions(self, partitions, threads=1, use_calculated_freqs=True, tree_search=True, batchsize=1, background=False):
        self.add_minsq_partitions(partitions)
        if isinstance(partitions, Partition):
            partitions = (partitions,)
        index_tuples = set(ix for partition in partitions for ix in partition.get_membership()).difference(
            self.lnl_cache.keys())
        if len(index_tuples) == 0:
            return

        # Collect argument list
        args = []
        to_delete = []
        try:
            for ix in index_tuples:
                conc = self.concatenate(ix)
                al = conc.alignment
                filename, delete = al.get_alignment_file(as_phylip=True)
                if delete:
                    to_delete.append(filename)
                partition = conc.qfile(default_dna="GTR", default_protein="LG", ml_freqs=True)
                tree = self.minsq_cache[ix]['tree']
                if use_calculated_freqs:
                    args.append((filename, partition, tree, tree_search, threads, PLL_RANDOM_SEED, conc.frequencies))
                else:
                    args.append((filename, partition, tree, tree_search, threads, PLL_RANDOM_SEED, None))

            # Distribute work
            with fileIO.TempFileList(to_delete):
                msg = 'Adding ML cluster trees'
                client = get_client()
                if client is None:
                    map_result = sequential_map(tasks.pll_task, args, msg)
                else:
                    map_result = parallel_map(client, tasks.pll_task, args, msg, batchsize, background)
                    if background:
                        return map_result

            # Process results
            pbar = setup_progressbar('Processing results', len(map_result))
            pbar.start()
            for i, (ix, result) in enumerate(zip(index_tuples, map_result)):
                self.lnl_cache[ix] = result
                pbar.update(i)
            pbar.finish()
        except:
            with fileIO.TempFileList(to_delete):
                pass

    def add_minsq_partitions(self, partitions, batchsize=1, background=False, **kwargs):
        if isinstance(partitions, Partition):
            partitions = (partitions,)
        index_tuples = set(ix for partition in partitions for ix in partition.get_membership()).difference(
            self.minsq_cache.keys())

        if len(index_tuples) == 0:
            return

        # Collect argument list
        args = []
        for ix in index_tuples:
            conc = self.concatenate(ix)
            args.append(conc.get_tree_collection_strings(**kwargs))

        # Distribute work
        msg = 'Adding MinSq cluster trees'
        client = get_client()
        if client is None:
            map_result = sequential_map(tasks.minsq_task, args, msg)
        else:
            map_result = parallel_map(client, tasks.minsq_task, args, msg, batchsize, background)
            if background:
                return map_result
        # Process results
        pbar = setup_progressbar('Processing results', len(map_result))
        pbar.start()
        for i, (ix, result) in enumerate(zip(index_tuples, map_result)):
            self.minsq_cache[ix] = result
            pbar.update(i)
        pbar.finish()

    def concatenate(self, index_tuple):
        """ Returns a Concatenation object that stitches together
        the alignments picked out by the index tuple """
        return self.collection.concatenate(index_tuple)

    def members(self, index_tuple):
        """ Gets records by their index, contained in the index_tuple """
        return [self.records[n] for n in index_tuple]

    def get_likelihood(self, partition, batchsize=1, **kwargs):
        """
        Return the sum of log-likelihoods for a partition.
        :param partition: Partition object
        :return: score (float)
        """
        indices = partition.get_membership()
        self.add_lnl_partitions(partition, batchsize=batchsize, **kwargs)
        results = [self.lnl_cache[ix] for ix in indices]
        return math.fsum(x['likelihood'] for x in results)

    def get_sse(self, partition, **kwargs):
        """
        Return the sum of squared errors score for a partition
        :param partition: Partition object
        :return: score (float)
        """
        indices = partition.get_membership()
        self.add_minsq_partitions(partition, **kwargs)
        results = [self.minsq_cache[ix] for ix in indices]
        return math.fsum(x['sse'] for x in results)

    # def get_fit(self, partition):
    #     """
    #                          harmonic mean of variances     sum of sq err
    #     Dimensionless fit ~  -------------------------- * ------------------
    #                            variance of distances      degrees of freedom
    #
    #     Return the dimensionless fit index for a partition
    #     in the tree.
    #     :param partition: Partition object
    #     :return: score (float)
    #     """
    #     results = self.get_results(partition, 'minsq')
    #     return math.fsum(x['fit'] for x in results)

    def simulate(self, partition, outdir, batchsize=1, **kwargs):
        """
        Simulate a set of alignments from the parameters inferred on a partition
        :param partition:
        :return:
        """
        indices = partition.get_membership()
        self.add_lnl_partitions(partition, **kwargs)
        results = [self.lnl_cache[ix] for ix in indices]
        places = dict((j,i) for (i,j) in enumerate(rec.name for rec in self.collection.records))

        # Collect argument list
        args = [None] * len(self.collection)
        for result in results:
            for partition in result['partitions'].values():
                place = places[partition['name']]
                args[place] = (len(self.collection[place]),
                               model_translate(partition['model']),
                               partition['frequencies'],
                               partition['alpha'],
                               result['ml_tree'],
                               partition['rates'] if 'rates' in partition else None)

        # Distribute work
        msg = 'Simulating'
        client = get_client()
        if client is None:
            map_result = sequential_map(client, tasks.simulate_task, args, msg)
        else:
            map_result = parallel_map(client, tasks.simulate_task, args, msg, batchsize, background)
            if background:
                return map_result

        # Process results
        for i, result in enumerate(map_result):
            orig = self.collection[i]
            simseqs = gapmask(result, orig.get_sequences())
            al = Alignment(simseqs, 'protein' if orig.is_protein() else 'dna')
            outfile = os.path.join(outdir, orig.name + '.phy')
            al.write_alignment(outfile, 'phylip', True)
