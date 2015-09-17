#!/usr/bin/env python
from __future__ import print_function

# standard lib
import glob
import hashlib
import itertools
import json
import math
import os
import random
import sys

# third party
import numpy as np
import phylo_utils
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
from .parutils import SequentialJobHandler
from .tree import Tree
from .utils import fileIO, setup_progressbar, model_translate, smooth_freqs
from .utils.decorators import lazyprop

import json
import sys
import logging
logger = logging.getLogger(__name__)

default_jobhandler = SequentialJobHandler()


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
    def __init__(self, file_format, input_dir):
        self.file_format = file_format
        self.input_dir = input_dir

    def __str__(self):
        msg = ('No records were found in {0} matching\n'
               '\tfile_format = {1}\n'.format(self.input_dir,
                                              self.file_format))
        return msg


class RecordsHandler(object):

    def __init__(
            self,
            records=None,
            input_dir=None,
            param_dir=None,
            file_format='phylip',
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
                                                header_grep)

        else:
            raise ValueError('Provide a list of records, '
                            'or the path to a set of alignments')

        if param_dir is not None:
            self.read_parameters(param_dir)

        if not self.records:
            raise NoRecordsError(file_format, input_dir)


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

    def read_alignments(self, input_dir, file_format, header_grep=None):
        """ Get list of alignment files from an input directory *.fa, *.fas and
        *.phy files only

        Stores in self.files """

        compression = ['', 'gz', 'bz2']

        if file_format == 'fasta':
            extensions = ['fa', 'fas', 'fasta']

        elif file_format == 'phylip':
            extensions = ['phy']

        else:
            extensions = []

        extensions = list('.'.join([x,y]) if y else x for x,y in itertools.product(extensions, compression))

        files = fileIO.glob_by_extensions(input_dir, extensions)
        files.sort(key=SORT_KEY)
        self._input_files = files
        records = []

        if self.show_progressbars:
            pbar = setup_progressbar("Loading files", len(files), simple_progress=True)
            pbar.start()

        for i, f in enumerate(files):
            if f.endswith('gz') or f.endswith('bz2'):
                with fileIO.TempFile() as tmpfile:
                    with fileIO.freader(f, f.endswith('gz'), f.endswith('bz2')) as reader, fileIO.fwriter(tmpfile) as writer:
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
            with fileIO.fwriter(os.path.join(output_dir, '{}.json'.format(rec.name)), gz=gz) as outfile:
                rec.parameters.write(outfile, indent=4)


class RecordsCalculatorMixin(object):

    def calc_distances(self, indices=None, task_interface=None, jobhandler=default_jobhandler, batchsize=1):
        """
        Calculate fast approximate intra-alignment pairwise distances and variances using
        ML (requires ML models to have been set up using `calc_trees`).
        :return: None (all side effects)
        """
        if indices is None:
            indices = list(range(len(self)))

        if task_interface is None:
            task_interface = tasks.MLDistanceTaskInterface()

        records = [self[i] for i in indices]

        # Assemble argument lists
        args, to_delete = task_interface.scrape_args(records)

        # Dispatch
        msg = '{} estimation'.format(task_interface.name)
        map_result = jobhandler(task_interface.get_task(), args, msg, batchsize)

        # Process results
        with fileIO.TempFileList(to_delete):
            # pbar = setup_progressbar('Processing results', len(map_result))
            # j = 0
            # pbar.start()
            for rec, result in zip(records, map_result):
                rec.parameters.partitions.distances = result['partitions'][0]['distances']
                rec.parameters.partitions.variances = result['partitions'][0]['variances']
                rec.parameters.nj_tree = result['nj_tree']
                # pbar.update(j+1)
                # j += 1
            # pbar.finish()

    def calc_trees(self, indices=None, task_interface=None, jobhandler=default_jobhandler, batchsize=1, **kwargs):
        if indices is None:
            indices = list(range(len(self)))

        if task_interface is None:
            task_interface = tasks.PllTaskInterface()

        records = [self[i] for i in indices]

        # Scrape args from records
        args, to_delete = task_interface.scrape_args(records, **kwargs)

        # Dispatch work
        msg = '{} Tree estimation'.format(task_interface.name)
        map_result = jobhandler(task_interface.get_task(), args, msg, batchsize)

        # Process results
        with fileIO.TempFileList(to_delete):
            for rec, result in zip(records, map_result):
                #logger.debug('Result - {}'.format(result))
                rec.parameters.construct_from_dict(result)

    def get_inter_tree_distances(self, metric, jobhandler=default_jobhandler, normalise=False, batchsize=1):
        """ Generate a distance matrix from a fully-populated Collection """
        metrics = {'euc': tasks.EuclideanTreeDistance,
                   'geo': tasks.GeodesicTreeDistance,
                   'rf': tasks.RobinsonFouldsTreeDistance,
                   'wrf': tasks.WeightedRobinsonFouldsTreeDistance}
        optioncheck(metric, metrics.keys())
        task_interface = metrics[metric]()
        args = task_interface.scrape_args(self.trees, normalise)
        msg = task_interface.name
        array = jobhandler(task_interface.get_task(), args, msg, batchsize)
        return DistanceMatrix.from_array(squareform(array), self.names)


class Collection(RecordsHandler, RecordsCalculatorMixin):
    """ Call:

    c = Collection(input_dir, file_format, datatype, tmpdir ...)
    c.calc_distances(), c.calc_TC_trees(), ...
    dm = c.distance_matrix('geo')
    cl = Clustering(dm)
    k = cl.spectral(4, prune='estimate', local_scale=7)
    p = Partition(k) """

    def num_species(self):
        """ Returns the number of species found over all records
        """
        all_headers = reduce(lambda x, y: set(x) | set(y),
                             (rec.get_names() for rec in self.records))
        return len(all_headers)

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

    def concatenate(self, indices):
        return Concatenation(self, indices)


class Scorer(object):

    def __init__(self, collection, cache_dir, task_interface):
        """
        Coordinates scoring of (usually) multilocus alignments in a partition
        """
        self.collection = collection
        self.cache_dir = cache_dir
        self.task_interface = task_interface
        self.cache = {}
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        if not os.path.exists(cache_dir):
            raise IOError('\'{}\' does not exist'.format(cache_dir))

    def get_id(self, grp):
        """
        Return a hash of the tuple of indices that specify the group
        """
        return self.cache.get(grp, hashlib.sha1(hex(hash(grp))).hexdigest())

    def check_work_done(self, grp):
        """
        Check for the existence of alignment and result files.
        """
        id_ = self.get_id(grp)
        concat_file = os.path.join(self.cache_dir, '{}.phy'.format(id_))
        result_file = os.path.join(self.cache_dir, '{}.{}.json'.format(id_, self.task_interface.name))
        return os.path.exists(concat_file), os.path.exists(result_file)

    def write_group(self, grp, overwrite=False, **kwargs):
        """
        Write the concatenated alignment to disk in the location specified by
        self.cache_dir
        """
        id_ = self.get_id(grp)
        alignment_done, result_done = self.check_work_done(grp)
        self.cache[grp] = id_
        al_filename = os.path.join(self.cache_dir, '{}.phy'.format(id_))
        qfile_filename = os.path.join(self.cache_dir, '{}.partitions.txt'.format(id_))
        if overwrite or not (alignment_done or result_done):
            conc = self.collection.concatenate(grp)
            al = conc.alignment
            al.write_alignment(al_filename, 'phylip', True)
            q = conc.qfile(**kwargs)
            with open(qfile_filename, 'w') as fl:
                fl.write(q + '\n')

    def get_group_result(self, grp, **kwargs):
        """
        Retrieve the results for a group. Needs this to already be calculated -
        errors out if result not available.
        """
        id_ = self.get_id(grp)
        self.cache[grp] = id_

        # Check if this file is already processed
        alignment_written, results_written = self.check_work_done(grp)

        if not results_written:
            if not alignment_written:
                self.write_group(grp, **kwargs)
            logger.error('Alignment {} has not been analysed - run analyse_cache_dir'.format(id_))
            raise ValueError('Missing result')
        else:
            with open(self.get_result_file(id_)) as fl:
                return json.load(fl)

    def get_result_file(self, id_):
        f = os.path.join(self.cache_dir, id_ + '.phy')
        return f.replace('.phy', '.{}.json'.format(self.task_interface.name))

    def write_partition(self, p, overwrite=False, **kwargs):
        for grp in p.get_membership():
            self.write_group(grp, overwrite, **kwargs)

    def analyse_cache_dir(self, jobhandler=None, batchsize=1, **kwargs):
        """
        Scan the cache directory and launch analysis for all unscored alignments
        using associated task handler. KWargs are passed to the tree calculating
        task managed by the TaskInterface in self.task_interface.
        Example kwargs:
          TreeCollectionTaskInterface: scale=1, guide_tree=None, 
                                       niters=10, keep_topology=False
          RaxmlTaskInterface: -------- partition_files=None, model=None, threads=1
          FastTreeTaskInterface: ----- No kwargs
        """
        if jobhandler is None:
            jobhandler = SequentialJobHandler()
        files = glob.glob(os.path.join(self.cache_dir, '*.phy'))
        #logger.debug('Files - {}'.format(files))
        records = []
        outfiles = []
        dna = self.collection[0].is_dna() # THIS IS ONLY A GUESS AT SEQ TYPE!!
        for infile in files:
            id_ = fileIO.strip_extensions(infile)
            outfile = self.get_result_file(id_)
            #logger.debug('Looking for {}: {}'.format(outfile, os.path.exists(outfile)))
            if not os.path.exists(outfile):
                record = Alignment(infile, 'phylip', True)
                records.append(record)
                outfiles.append(outfile)

        if len(records) == 0:
            return []

        args, to_delete = self.task_interface.scrape_args(records, outfiles=outfiles, **kwargs)
        logger.debug('Args - {}'.format(args))

        with fileIO.TempFileList(to_delete):
            result = jobhandler(self.task_interface.get_task(), args, 'Cache dir analysis', batchsize)
            for (out, res) in zip(outfiles, result):
                if not os.path.exists(out) and res:
                    with open(out, 'w') as outfl:
                        json.dump(res, outfl)

        return result

    def get_partition_score(self, p):
        """
        Assumes analysis is done and written to id.json!
        """
        scores = []
        for grp in p.get_membership():
            try:
                result = self.get_group_result(grp)
                scores.append(result['likelihood'])
            except ValueError:
                scores.append(None)
        return sum(scores)

    def get_partition_trees(self, p):
        """
        Return the trees associated with a partition, p
        """
        trees = []
        for grp in p.get_membership():
            try:
                result = self.get_group_result(grp)
                trees.append(result['ml_tree'])
            except ValueError:
                trees.append(None)
                logger.error('No tree found for group {}'.format(grp))
        return trees

    def get_partition_members(self, p):
        result = []
        for grp in p.get_membership():
            members = [self.collection[i].name for i in grp]
            result.append(members)
        return result

    def get_partition_results(self, p):
        results = []
        for grp in p.get_membership():
            try:
                result = self.get_group_result(grp)
                results.append(result)
            except ValueError:
                results.append(None)
                logger.error('No result found for group {}'.format(grp))
        return results

    def clean_cache(self):
        files = glob.glob(os.path.join(self.cache_dir, '*.json'))
        for f in files:
            id_ = fileIO.strip_extensions(f)
            alfile = os.path.join(self.cache_dir, '{}.phy'.format(id_))
            qfile = os.path.join(self.cache_dir, '{}.partitions.txt'.format(id_))
            if os.path.exists(alfile): os.remove(alfile)
            if os.path.exists(qfile): os.remove(qfile)
            alfile = os.path.join(self.cache_dir, '{}.phy.reduced'.format(id_))
            qfile = os.path.join(self.cache_dir, '{}.partitions.txt.reduced'.format(id_))
            if os.path.exists(alfile): os.remove(alfile)
            if os.path.exists(qfile): os.remove(qfile)

    def simulate(self, partition, outdir, jobhandler=default_jobhandler, batchsize=1, **kwargs):
        """
        Simulate a set of alignments from the parameters inferred on a partition
        :param partition:
        :return:
        """

        results = self.get_partition_results(partition)
        DEFAULT_DNA_MODEL = 'GTR'
        DEFAULT_PROTEIN_MODEL = 'LG08'

        # Collect argument list
        args = [None] * len(self.collection)

        for result in results:
            if len(result['partitions']) > 1:
                places = dict((j,i) for (i,j) in enumerate(rec.name for rec in self.collection.records))
                for partition in result['partitions'].values():
                    place = places[partition['name']]
                    model = partition.get('model')
                    freqs = partition.get('frequencies')
                    rates = partition.get('rates')
                    alpha = partition.get('alpha')
                    tree  = str(result['ml_tree'])
                    if model is None:
                        model = DEFAULT_DNA_MODEL if self.collection[place].is_dna() else DEFAULT_PROTEIN_MODEL
                    if freqs is not None:
                        freqs = smooth_freqs(freqs)
                    args[place] = (len(self.collection[place]),
                                   model_translate(model),
                                   freqs,
                                   alpha,
                                   tree,
                                   rates)
            else:
                model = result['partitions']['0'].get('model')
                freqs = result['partitions']['0'].get('frequencies')
                rates = result['partitions']['0'].get('rates')
                alpha = result['partitions']['0'].get('alpha')
                tree  = str(result['ml_tree'])
                if freqs is not None:
                    freqs = smooth_freqs(freqs)
                use_default_model = (model is None)
                for i in range(len(self.collection)):
                    if use_default_model:
                        model = DEFAULT_DNA_MODEL if self.collection[i].is_dna() else DEFAULT_PROTEIN_MODEL
                    args[i] = (len(self.collection[i]),
                               model_translate(model),
                               freqs,
                               alpha,
                               tree,
                               rates)

        # Distribute work
        msg = 'Simulating'
        map_result = jobhandler(tasks.simulate_task, args, msg, batchsize)

        # Process results
        for i, result in enumerate(map_result):
            orig = self.collection[i]
            simseqs = gapmask(result, orig.get_sequences())
            al = Alignment(simseqs, 'protein' if orig.is_protein() else 'dna')
            outfile = os.path.join(outdir, orig.name + '.phy')
            al.write_alignment(outfile, 'phylip', True)

class Optimiser(object):
    def __init__(self, scorer):
        self.scorer = scorer
        self.insts = init_perlocus_likelihood_objects()

    def init_perlocus_likelihood_objects(self):
        from pllpy import pll
        c = self.scorer.collection
        insts = []
        to_delete = []
        for rec in c:
            alignment_file, delete = rec.get_alignment_file()
            if delete:
                to_delete.append(alignment_file)
            model = 'GTR' if rec.is_dna() else 'LG'
            partition_string = '{}, {} = 1 - {}'.format(model, rec.name, len(rec))
            inst = pll(alignment_file, partition_string, rec.tree, 1, 12345)
            insts.append(inst)
        return insts
