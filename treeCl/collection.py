#!/usr/bin/env python
from __future__ import print_function

# standard lib
import json
import math
import os
import sys
import random
import time
import timeit

# third party

# treeCl
from treeCl.concatenation import Concatenation
from treeCl.clustering import Partition
from treeCl.tasks import tasks
from treeCl.tasks.celery import app
from distance_matrix import DistanceMatrix
from alignment import Alignment
from parameters import PartitionParameters
from utils import fileIO, setup_progressbar
from errors import optioncheck, directorycheck
from constants import SORT_KEY, PLL_RANDOM_SEED

DISTRIBUTED_TASK_QUEUE_INSPECT = app.control.inspect()


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
    ):

        self._records = None
        self._input_files = None

        if records is not None:
            self.records = records

        elif input_dir is not None:
            input_dir = os.path.abspath(input_dir)
            directorycheck(input_dir)
            optioncheck(file_format, ['fasta', 'phylip'])
            self.records = self.read_alignments(input_dir,
                                                file_format,
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

    @property
    def trees(self):
        """ Returns a list of trees in SORT_KEY order """
        try:
            return [rec.tree for rec in self]
        except ValueError:
            return []

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

            record.name = (fileIO.strip_extensions(f))
            records.append(record)
            pbar.update(i)
        pbar.finish()
        return records

    def read_parameters(self, input_dir):
        """ Read a directory full of tree files, matching them up to the
        already loaded alignments """

        pbar = setup_progressbar("Loading parameters", len(self.records))
        pbar.start()
        for i, rec in enumerate(self.records):
            try:
                with open(os.path.join(input_dir, '{}.json'.format(rec.name))) as infile:
                    parameters = json.load(infile, parse_int=True)
                    rec.parameters = parameters
                    if 'partitions' in rec.parameters:
                        rec.parameters['partitions'] = {int(k): v for (k, v) in rec.parameters['partitions'].iteritems()}
                        try:
                            freqs = rec.parameters['partitions'][0]['frequencies']
                            alpha = rec.parameters['partitions'][0]['alpha']
                            rec.set_substitution_model('GTR' if rec.is_dna() else 'LG08')
                            rec.set_gamma_rate_model(4, alpha)
                            rec.set_frequencies(freqs)
                            if rec.is_dna():
                                rec.set_rates(result['partitions'][0]['rates'], 'ACGT')
                        except KeyError:
                            pass

                        try:
                            dists = rec.parameters['partitions'][0]['distances']
                            rec.set_distance_matrix(dists)
                        except KeyError:
                            pass

            except IOError:
                continue
            finally:
                pbar.update(i)
        pbar.finish()

    def write_parameters(self, output_dir):
        if not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
            except IOError as err:
                sys.stderr.write(err.message)
                raise err

        for rec in self.records:
            with fileIO.fwriter(os.path.join(output_dir, '{}.json'.format(rec.name)), gz=True) as outfile:
                rec.parameters.write(outfile, indent=4)


####### TASKS ##########################################################################################################

    def fast_calc_distances(self):
        if DISTRIBUTED_TASK_QUEUE_INSPECT.active() is None:
            self._fast_calc_distances_sequential()
        else:
            self._fast_calc_distances_async()

    # noinspection PyUnresolvedReferences
    def _fast_calc_distances_sequential(self):
        """ Calculates within-alignment pairwise distances and variances for every
        alignment. Uses fast Jukes-Cantor method.
        :return: void"""
        pbar = setup_progressbar('Calculating fast approximate distances', len(self))
        pbar.start()
        for i, rec in enumerate(self.records):
            rec.fast_compute_distances()
            pbar.update(i)
            params = PartitionParameters()
            params.distances = rec.get_distances().tolist()
            params.variances = rec.get_variances().tolist()
            rec.parameters.partitions.append(params)
            rec.parameters.nj_tree = rec.get_bionj_tree()
        pbar.finish()

    # noinspection PyUnresolvedReferences
    def _fast_calc_distances_async(self):
        from celery import group
        jobs = []
        to_delete = []
        for rec in self:
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            jobs.append((filename))

        with fileIO.TempFileList(to_delete):
            job_group = group(tasks.fast_calc_distances_task.s(args) for args in jobs)()
            pbar = setup_progressbar('Calculating fast distances', len(jobs), simple_progress=True)
            pbar.start()
            while not job_group.ready():
                time.sleep(2)
                pbar.update(job_group.completed_count())
            pbar.finish()

        pbar = setup_progressbar('Processing results', len(jobs))
        j = 0
        pbar.start()
        for i, async_result in enumerate(job_group.results):
            rec = self[i]
            result = async_result.get()
            distances = result['distances']
            variances = result['variances']
            tree = result['tree']
            rec.parameters.nj_tree = tree
            params = PartitionParameters()
            params.distances = distances
            params.variances = variances
            rec.parameters.partitions = [params]
            pbar.update(i)
        pbar.finish()

    def calc_distances(self):
        """ Calculates within-alignment pairwise distances and variances for every
        alignment. Uses slow ML optimisation, and depends on the Alignment record
        having had appropriate ML parametric models set up in advance.
        :return: void
        """
        if DISTRIBUTED_TASK_QUEUE_INSPECT.active() is None:
            self._calc_distances_sequential()
        else:
            self._calc_distances_async()

    def _calc_distances_sequential(self):
        pbar = setup_progressbar('Calculating ML distances', len(self))
        pbar.start()
        to_delete = []
        for i, rec in enumerate(self.records):
            # Get file
            filename, delete = rec.get_alignment_file()
            if delete:
                to_delete.append(filename)
            # Get input dict
            model = {'partitions': {}}
            data = {'alpha': rec.parameters.partitions.alpha, 'frequencies': rec.parameters.partitions.frequencies}
            if rec.is_dna():
                data['rates'] = rec.parameters.partitions.rates
            model['partitions'][0] = data
            # Launch local task
            result = tasks.calc_distances_task(model, filename)
            rec.parameters.partitions.distances = result['partitions'][0]['distances']
            rec.parameters.partitions.variances = result['partitions'][0]['variances']
            rec.parameters.nj_tree = result['nj_tree']
            pbar.update(i)
        with fileIO.TempFileList(to_delete):
            pbar.finish()

    def _calc_distances_async(self):
        from celery import group

        jobs = []
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
            jobs.append((model, filename))

        with fileIO.TempFileList(to_delete):
            job_group = group(tasks.calc_distances_task.subtask(args) for args in jobs)()
            pbar = setup_progressbar('Calculating ML distances', len(jobs))
            pbar.start()
            while not job_group.ready():
                time.sleep(2)
                pbar.update(job_group.completed_count())
            pbar.finish()

        pbar = setup_progressbar('Processing results', len(jobs))
        j = 0
        pbar.start()
        for i, async_result in enumerate(job_group.results):
            result = async_result.get(timeout=20)
            rec = self[i]
            rec.parameters.partitions.distances = result['partitions'][0]['distances']
            rec.parameters.partitions.variances = result['partitions'][0]['variances']
            rec.parameters.nj_tree = result['nj_tree']
            pbar.update(j+1)
            j += 1

    def calc_trees(self, threads=1, indices=None):
        if DISTRIBUTED_TASK_QUEUE_INSPECT.active() is None:
            self._calc_trees_sequential(threads, indices)
        else:
            self._calc_trees_async(threads, indices)

    def _calc_trees_sequential(self, threads=1, indices=None):
        """ Use pllpy to calculate maximum-likelihood trees
        :return: void
        """

        if indices is None:
            indices = list(range(len(self)))
        else:
            indices = indices

        pbar = setup_progressbar('Calculating ML trees', len(indices))
        pbar.start()

        to_delete = []
        for i, rec in enumerate(self[i] for i in indices):
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            partition = '{}, {} = 1 - {}'.format('DNA' if rec.is_dna() else 'LGX', rec.name, len(rec))
            try:
                tree = rec.tree
            except AttributeError:
                tree = True
            result = tasks.pll_task(filename, partition, tree, threads, PLL_RANDOM_SEED)
            rec.set_params_from_pll_result(result)
            pbar.update(i)

        with fileIO.TempFileList(to_delete):
            pbar.finish()

    # noinspection PyUnresolvedReferences
    def _calc_trees_async(self, threads=1, indices=None, allow_retry=True):
        """ Use pllpy to calculate maximum-likelihood trees, and use celery to distribute
        the computation across cores
        :return: void
        """
        from celery import group
        from celery.exceptions import TimeoutError
        if indices is None:
            indices = list(range(len(self)))
        jobs = []
        to_delete = []
        for i in indices:
            rec = self[i]
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            partition = '{}, {} = 1 - {}'.format('DNA' if rec.is_dna() else 'LGX', rec.name, len(rec))
            tree = rec.parameters.nj_tree if rec.parameters.nj_tree is not None else True
            jobs.append((filename, partition, tree, threads, PLL_RANDOM_SEED))

        with fileIO.TempFileList(to_delete):
            job_group = group(tasks.pll_task.subtask(args) for args in jobs)()
            pbar = setup_progressbar('Calculating ML trees', len(jobs))
            pbar.start()
            while not job_group.ready():
                time.sleep(2)
                pbar.update(job_group.completed_count())
            pbar.finish()

        pbar = setup_progressbar('Processing results', len(jobs))
        j = 0
        pbar.start()
        retries = []
        for i, async_result in zip(indices, job_group.results):
            try:
                result = async_result.get(timeout=20)
            except TimeoutError:
                retries.append(i)
            rec = self[i]
            rec.set_params_from_pll_result(result)
            pbar.update(j+1)
            j += 1
        if retries > [] and allow_retry:
            self._calc_trees_async(1, retries, False)

    def get_inter_tree_distances(self, metric, **kwargs):
        """ Generate a distance matrix from a fully-populated Collection """
        distribute_tasks = DISTRIBUTED_TASK_QUEUE_INSPECT.active() is not None
        return DistanceMatrix(self.trees, metric, distribute_tasks=distribute_tasks, **kwargs)

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

    def add_lnl_partitions(self, partitions, threads=1):
        self.add_minsq_partitions(partitions)
        if isinstance(partitions, Partition):
            partitions = (partitions,)
        index_tuples = set(ix for partition in partitions for ix in partition.get_membership()).difference(
            self.lnl_cache.keys())
        if len(index_tuples) > 0:
            if DISTRIBUTED_TASK_QUEUE_INSPECT.active is None:
                self._add_lnl_sequential(index_tuples, threads)
            else:
                self._add_lnl_async(index_tuples, threads)

    def _add_lnl_sequential(self, index_tuples, threads=1):
        pbar = setup_progressbar('Adding ML cluster trees: ', len(index_tuples))
        pbar.start()

        to_delete = []
        for i, ix in enumerate(index_tuples):
            conc = self.concatenate(ix)
            al = conc.alignment
            filename, delete = al.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            partition = conc.qfile(dna_model="GTR", protein_model="LG", ml_freqs=True)
            tree = self.minsq_cache[ix]['tree']
            self.lnl_cache[ix] = tasks.pll_task(filename, partition, tree, threads, PLL_RANDOM_SEED)
            pbar.update(i)

        with fileIO.TempFileList(to_delete):
            pbar.finish()

    def _add_lnl_async(self, index_tuples, threads=1):
        from celery import group

        jobs = []
        to_delete = []
        for ix in index_tuples:
            conc = self.concatenate(ix)
            al = conc.alignment
            filename, delete = al.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            partition = conc.qfile(dna_model="GTR", protein_model="LG", ml_freqs=True)
            tree = self.minsq_cache[ix]['tree']
            jobs.append((filename, partition, tree, threads, PLL_RANDOM_SEED))

        with fileIO.TempFileList(to_delete):
            job_group = group(tasks.pll_task.subtask(args) for args in jobs)()
            pbar = setup_progressbar('Adding ML cluster trees', len(jobs))
            pbar.start()
            while not job_group.ready():
                time.sleep(2)
                pbar.update(job_group.completed_count())
            pbar.finish()

        pbar = setup_progressbar('Processing results', len(jobs))
        pbar.start()
        for i, (ix, async_result) in enumerate(zip(index_tuples, job_group.results)):
            self.lnl_cache[ix] = async_result.get(timeout=20)
            pbar.update(i)
        pbar.finish()

    def add_minsq_partitions(self, partitions):
        if isinstance(partitions, Partition):
            partitions = (partitions,)
        index_tuples = set(ix for partition in partitions for ix in partition.get_membership()).difference(
            self.lnl_cache.keys())
        if len(index_tuples) > 0:
            if DISTRIBUTED_TASK_QUEUE_INSPECT.active is None:
                self._add_minsq_sequential(index_tuples)
            else:
                self._add_minsq_async(index_tuples)

    def _add_minsq_sequential(self, index_tuples):
        pbar = setup_progressbar('Adding MinSq cluster trees: ', len(index_tuples))
        pbar.start()
        for i, ix in enumerate(index_tuples):
            conc = self.concatenate(ix)
            dv, lab, gm, tree = conc.get_tree_collection_strings()
            self.minsq_cache[ix] = tasks.minsq_task(dv, lab, gm, tree)
            pbar.update(i)
        pbar.finish()

    def _add_minsq_async(self, index_tuples):
        from celery import group

        jobs = []
        for ix in index_tuples:
            conc = self.concatenate(ix)
            jobs.append(conc.get_tree_collection_strings())

        job_group = group(tasks.minsq_task.subtask(args) for args in jobs)()
        pbar = setup_progressbar('Adding MinSq cluster trees', len(jobs))
        pbar.start()
        while not job_group.ready():
            time.sleep(2)
            pbar.update(job_group.completed_count())
        pbar.finish()

        pbar = setup_progressbar('Processing results', len(jobs))
        pbar.start()
        for i, (ix, async_result) in enumerate(zip(index_tuples, job_group.results)):
            self.minsq_cache[ix] = async_result.get(timeout=20)
            pbar.update(i)
        pbar.finish()

    def get_lnl_partition(self, partition):
        """ Calculates concatenated trees for a Partition """
        index_tuple_list = partition.get_membership()

        return [self._get_lnl(ix, nthreads) for ix in index_tuple_list]


    def _get_lnl(self, index_tuple):
        """
        Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree
        :param index_tuple: tuple of indexes at which to find the alignments.
        Result is memoized in self.lnl_cache dictionary
        :return: dictionary
        """
        try:
            return self.lnl_cache[index_tuple]
        except KeyError:
            conc = self.concatenate(index_tuple)
            partitions = conc.qfile(dna_model="GTR", protein_model="LGX")
            tree = self._get_minsq(index_tuple)['tree']
            result = self.lnl_cache[index_tuple] = \
                tasks.pll_task()#todo
            return result

    def get_minsq_partition(self, partition):
        """ Calculates concatenated trees for a Partition """
        index_tuple_list = partition.get_membership()
        return [self._get_minsq(index_tuple) for index_tuple in index_tuple_list]

    def _get_minsq_partition_sequential(self):
        pass

    def _get_minsq_partition_async(self):
        pass

    def _get_minsq(self, index_tuple):
        """
        Takes a tuple of indices. Concatenates the records in the record
        list at these indices, and builds a tree. Returns the tree
        Result is memoized in self.minsq_cache dictionary
        :param index_tuple: tuple of indexes at which to find the alignments.
        :return: dictionary
        """
        try:
            return self.minsq_cache[index_tuple]
        except KeyError:
            conc = self.concatenate(index_tuple)
            tree, sse = conc.minsq_tree()
            tree.deroot()
            n_tips = len(tree)
            result = dict(tree=tree.newick, sse=sse, fit=sse / (2 * (n_tips - 2) * (n_tips - 3)), names=conc.names)
            self.minsq_cache[index_tuple] = result
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

    def get_results(self, partition, criterion, use_celery=False, nthreads=1):
        """
        Return the results for scoring a partition - either the sum of log likelihoods,
        or the total min squares dimensionless fit index
        :param partition: Partition object
        :param criterion: either 'minsq' or 'lnl'
        :return: score (float)
        """
        optioncheck(criterion, ['lnl', 'minsq'])
        results = (self.get_lnl_partition(partition, use_celery, nthreads) if criterion == 'lnl'
                   else self.get_minsq_partition(partition))
        return results

    def get_likelihood(self, partition, use_celery=False, nthreads=1):
        """
        Return the sum of log-likelihoods for a partition.
        :param partition: Partition object
        :return: score (float)

        """
        results = self.get_results(partition, 'lnl', use_celery, nthreads)
        return math.fsum(x['likelihood'] for x in results)

    def get_sse(self, partition):
        """
        Return the sum of squared errors score for a partition
        :param partition: Partition object
        :return: score (float)
        """
        results = self.get_results(partition, 'minsq')
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
