#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from builtins import next
from builtins import hex
from builtins import str
from builtins import zip
from builtins import range
from builtins import object

# standard lib
import glob
import hashlib
import itertools
import json
import os
import random
import sys
import tempfile
from functools import reduce

# third party
import numpy as np
import phylo_utils
from scipy.spatial.distance import squareform
from tree_distance import PhyloTree

# treeCl
from .alignment import Alignment
from .concatenation import Concatenation
from .constants import SORT_KEY, ISPY3
from .distance_matrix import DistanceMatrix
from .errors import optioncheck, directorycheck
from . import tasks
from .partition import Partition
from .parutils import SequentialJobHandler
from .utils import fileIO, setup_progressbar, model_translate, smooth_freqs, create_gamma_model, flatten_list
from .utils.decorators import lazyprop

# set up logging
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


def transform_fn(table, amount=2.0):
    tmp = table**(1.0/amount)
    return tmp / tmp.sum(1)[:,np.newaxis]


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
            trees_dir=None,
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
            self.input_dir = input_dir

        else:
            raise ValueError('Provide a list of records, '
                            'or the path to a set of alignments')

        if param_dir is not None and trees_dir is None:
            self.read_parameters(param_dir)

        elif trees_dir is not None and param_dir is None:
            self.read_trees(trees_dir)

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
            if f.endswith('.gz') or f.endswith('.bz2'):
                fd, tmpfile = tempfile.mkstemp()
                with fileIO.freader(f, f.endswith('.gz'), f.endswith('.bz2')) as reader,\
                     fileIO.fwriter(tmpfile) as writer:
                    for line in reader:
                        if ISPY3:
                            line = line.decode()
                        writer.write(line)
                try:
                    record = Alignment(tmpfile, file_format, True)
                except ValueError:
                    record = Alignment(tmpfile, file_format, False)
                finally:
                    os.close(fd)
                    os.unlink(tmpfile)

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

    def read_trees(self, input_dir):
        """ Read a directory full of tree files, matching them up to the
        already loaded alignments """

        if self.show_progressbars:
            pbar = setup_progressbar("Loading trees", len(self.records))
            pbar.start()

        for i, rec in enumerate(self.records):
            hook = os.path.join(input_dir, '{}.nwk*'.format(rec.name))
            filename = glob.glob(hook)
            try:
                with fileIO.freader(filename[0]) as infile:
                    tree = infile.read().decode('utf-8')

                d = dict(ml_tree=tree)

                rec.parameters.construct_from_dict(d)

            except (IOError, IndexError):
                continue

            finally:
                if self.show_progressbars:
                    pbar.update(i)

        if self.show_progressbars:
            pbar.finish()

    def read_parameters(self, input_dir):
        """ Read a directory full of json parameter files, matching them up to the
        already loaded alignments """

        if self.show_progressbars:
            pbar = setup_progressbar("Loading parameters", len(self.records))
            pbar.start()
        for i, rec in enumerate(self.records):
            hook = os.path.join(input_dir, '{}.json*'.format(rec.name))
            filename = glob.glob(hook)
            try:
                with fileIO.freader(filename[0]) as infile:
                    d = json.loads(infile.read().decode('utf-8'), parse_int=True)

                rec.parameters.construct_from_dict(d)

            except (IOError, IndexError):
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

    def calc_trees(self, indices=None, task_interface=None, jobhandler=default_jobhandler, batchsize=1,
                   **kwargs):
        """
        Infer phylogenetic trees for the loaded Alignments

        :param indices: Only run inference on the alignments at these given indices
        :param task_interface: Inference tool specified via TaskInterface (default RaxmlTaskInterface)
        :param jobhandler: Launch jobs via this JobHandler (default SequentialJobHandler; also available are
            ThreadpoolJobHandler and ProcesspoolJobHandler for running inference in parallel)
        :param batchsize: Batch size for Thread- or ProcesspoolJobHandlers)
        :param kwargs: Remaining arguments to pass to the TaskInterface
        :return: None
        """
        if indices is None:
            indices = list(range(len(self)))

        if task_interface is None:
            task_interface = tasks.RaxmlTaskInterface()

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

    def get_inter_tree_distances(self, metric, jobhandler=default_jobhandler, normalise=False, min_overlap=4, batchsize=1):
        """ Generate a distance matrix from a fully-populated Collection """
        metrics = {'euc': tasks.EuclideanTreeDistance,
                   'geo': tasks.GeodesicTreeDistance,
                   'rf': tasks.RobinsonFouldsTreeDistance,
                   'wrf': tasks.WeightedRobinsonFouldsTreeDistance,
                   'fasteuc': tasks.EqualLeafSetEuclideanTreeDistance,
                   'fastgeo': tasks.EqualLeafSetGeodesicTreeDistance,
                   'fastrf': tasks.EqualLeafSetRobinsonFouldsTreeDistance,
                   'fastwrf': tasks.EqualLeafSetWeightedRobinsonFouldsTreeDistance}
        optioncheck(metric, list(metrics.keys()))
        task_interface = metrics[metric]()
        if metric.startswith('fast'):
            trees = (PhyloTree(newick, False) for newick in self.trees)
        else:
            trees = self.trees
        args = task_interface.scrape_args(trees, normalise, min_overlap)
        logger.debug('{}'.format(args))
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

    def species_set(self):
        return reduce(lambda x, y: set(x) | set(y),
                      (rec.get_names() for rec in self.records))

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
            return [next(iterable) for _ in range(n)]

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
        thehash = hex(hash(grp))
        if ISPY3:  # use default encoding to get bytes
            thehash = thehash.encode()
        return self.cache.get(grp, hashlib.sha1(thehash).hexdigest())

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
        # logger.debug('Args - {}'.format(args))

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
            simseqs = gapmask(result.items(), orig.get_sequences())
            al = Alignment(simseqs, alphabet=('protein' if orig.is_protein() else 'dna'))
            outfile = os.path.join(outdir, orig.name + '.phy')
            al.write_alignment(outfile, 'phylip', True)


class Optimiser(object):
    """ Perform the Classification-Expectation-Maximisation (CEM) algorithm (1) 
        Brief:
        Optimise the assignment of N data points to K groups by cycling through 3 steps
        - Expectation    - Calculate probabilities of group membership for each data point x_i
                           to group P_k, based on current model parameters
        - Classification - Assign each data point x_i to one group P_k according to the
                           probability of membership
        - Maximisation   - Update model parameters according to new membership

        (1) Celeux,G. and Govaert,G. (1992) A classification EM algorithm for clustering 
        and two stochastic versions. Comput. Stat. Data Anal.,14,315-332"""

    def __init__(self, scorer, numgrp, partition=None, **kwargs):
        self.scorer = scorer
        self.numgrp = numgrp
        self.partition = None
        self.prev_partition = None
        if partition is not None:
            self.set_partition(partition)
        self.insts = self.init_perlocus_likelihood_objects(**kwargs)
        self.names_to_indices = dict((rec.name, i) for (i, rec) in enumerate(scorer.collection))
        self.iterations = 0
        self.log = []
        self.lktable = None
        self.table = None

    def expect(self, use_proportions=True):
        """ The Expectation step of the CEM algorithm """
        changed = self.get_changed(self.partition, self.prev_partition)
        lk_table = self.generate_lktable(self.partition, changed, use_proportions)
        self.table = self.likelihood_table_to_probs(lk_table)

    def classify(self, table, weighted_choice=False, transform=None):
        """ The Classification step of the CEM algorithm """
        assert table.shape[1] == self.numgrp
        if weighted_choice:
            if transform is not None:
                probs = transform_fn(table.copy(), transform)  #
            else:
                probs = table.copy()
            cmprobs = probs.cumsum(1)
            logger.info('Probabilities\n{}'.format(probs))
            r = np.random.random(cmprobs.shape[0])
            search = np.apply_along_axis(np.searchsorted, 1, cmprobs, r) # Not very efficient
            assignment = np.diag(search)
        else:
            probs = table
            assignment = np.where(probs==probs.max(1)[:, np.newaxis])[1]
        logger.info('Assignment\n{}'.format(assignment))
        assignment = self._fill_empty_groups(probs, assignment)  # don't want empty groups
        new_partition = Partition(tuple(assignment))
        self.set_partition(new_partition)

    def maximise(self, **kwargs):
        """ The Maximisation step of the CEM algorithm """
        self.scorer.write_partition(self.partition)
        self.scorer.analyse_cache_dir(**kwargs)
        self.likelihood = self.scorer.get_partition_score(self.partition)
        self.scorer.clean_cache()
        changed = self.get_changed(self.partition, self.prev_partition)
        self.update_perlocus_likelihood_objects(self.partition, changed)
        return self.partition, self.likelihood, sum(inst.get_likelihood() for inst in self.insts)

    def iterate(self, use_proportions=True, weighted_choice=False, transform=None, **kwargs):
        self.expect(use_proportions)
        self.classify(self.table, weighted_choice, transform)
        self.iterations += 1
        result = self.maximise(**kwargs)
        self.log.append(result)
        return result

    def random_partition(self, ngroups):
        items = len(self.scorer.collection)
        r = np.zeros(items)
        r[:ngroups] = np.arange(ngroups)
        r[ngroups:] = np.random.randint(ngroups, size=items-ngroups)
        np.random.shuffle(r)
        return Partition(tuple(r))

    def set_partition(self, partition):
        """
        Store the partition in self.partition, and
        move the old self.partition into self.prev_partition
        """
        assert len(partition) == self.numgrp
        self.partition, self.prev_partition = partition, self.partition

    def init_perlocus_likelihood_objects(self, **kwargs):
        c = self.scorer.collection
        insts = []
        for rec in c:
            gamma = create_gamma_model(rec, list(c.species_set() - set(rec.get_names())), **kwargs)
            gamma.set_tree(rec.tree)
            insts.append(gamma)
        return insts

    def get_cluster_at_index(self, i):
        """ 
        Return the cluster membership of locus i, according to current
        assignment 
        """
        return self.partition.partition_vector[i]

    def get_changed(self, p1, p2):
        """ 
        Return the loci that are in clusters that have changed between
        partitions p1 and p2 
        """
        if p1 is None or p2 is None:
            return list(range(len(self.insts)))
        return set(flatten_list(set(p1) - set(p2)))

    def update_perlocus_likelihood_objects_old(self, partition, changed):
        results = self.scorer.get_partition_results(partition)
        UNPARTITIONED=False
        for result in results:
            subdict = result['partitions']
            if len(subdict) > 1:
                for k in subdict:
                    p = subdict[k]
                    index = self.names_to_indices[p['name']]
                    if index in changed:
                        inst = self.insts[index]
                        self._update_likelihood_model(inst, p, result['ml_tree'])
            else:
                UNPARTITIONED=True # not nice, but I'm in a hurry

        if UNPARTITIONED:
            prev_partition = self.partition

            for i in changed:
                cluster = self.get_cluster_at_index(i)
                inst = self.insts[i]
                result = results[cluster]
                p = result['partitions']['0']
                self._update_likelihood_model(inst, p, result['ml_tree'])

    def update_perlocus_likelihood_objects(self, partition, changed):
        for grp in partition:
            result = self.scorer.get_group_result(grp)
            tree = result['ml_tree']
            for i in grp:
                self._update_likelihood_model(self.insts[i], result['partitions']['0'], tree)
                # self.insts[i].set_tree(tree)
                # self.insts[i].update_alpha(result['partitions']['0']['alpha'])

    def _update_likelihood_model(self, inst, partition_parameters, tree):
        """ 
        Set parameters of likelihood model - inst -
        using values in dictionary - partition_parameters -,
        and - tree -
        """
        # Build transition matrix from dict
        model = partition_parameters['model']
        freqs = partition_parameters.get('frequencies')
        if model == 'LG':
            subs_model = phylo_utils.models.LG(freqs)
        elif model == 'WAG':
            subs_model = phylo_utils.models.WAG(freqs)
        elif model == 'GTR':
            rates = partition_parameters.get('rates')
            subs_model = phylo_utils.models.GTR(rates, freqs, True)
        else:
            raise ValueError("Can't handle this model: {}".format(model))
        tm = phylo_utils.markov.TransitionMatrix(subs_model)
        
        # Read alpha value
        alpha = partition_parameters['alpha']
        inst.set_tree(tree)
        inst.update_alpha(alpha)
        inst.update_transition_matrix(tm)

    def generate_lktable(self, partition, changed, use_proportions=True):
        trees = self.scorer.get_partition_trees(partition)

        # Try to call up table from previous step
        prev_lktable = self.lktable

        if use_proportions:
            sizes = np.array([len(x) for x in partition.get_membership()], dtype=np.double)
            total = partition.num_elements()
            logproportions = np.log(sizes/total)
        else:
            logproportions = np.zeros(partition.num_elements())

        lktable = np.zeros((len(self.insts), self.numgrp))
        for i, gamma in enumerate(self.insts):
            if i in changed or prev_lktable is None:
                for j, t in enumerate(trees):
                    gamma.set_tree(t)
                    lktable[i, j] = logproportions[j] + gamma.get_likelihood()
            else:
                lktable[i] = prev_lktable[i]
        self.lktable = lktable
        return lktable

    def likelihood_table_to_probs(self, lktable):
        """
        Calculates this formula (1), given the log of the numerator as input
                     
                     p_k * f(x_i, a_k)
        t_k(x_i) = -----------------------
                    ---K
                    \   p_k * f(x_i, a_k)
                    /__k=1
        
        x_i is data point i
        P_k is cluster k of K
        t_k is the posterior probability of x_i belonging to P_k
        p_k is the prior probability of belong to P_k (the proportional size of P_k)
        f(x, a) is the likelihood of x with parameters a
        """
        m = lktable.max(1)  # row max of lktable
        shifted = lktable-m[:,np.newaxis]  # shift lktable of log-likelihoods to a non-underflowing range
        expsum = np.exp(shifted).sum(1)  # convert logs to (scaled) normal space, and sum the rows
        logexpsum = np.log(expsum)+m  # convert back to log space, and undo the scaling
        return np.exp(lktable - logexpsum[:, np.newaxis])

    def _fill_empty_groups(self, probs, assignment):
        new_assignment = np.array(assignment.tolist())
        for k in range(probs.shape[1]):
            if np.count_nonzero(assignment==k) == 0:
                logger.info('Group {} became empty'.format(k))
                # Group k is empty, so needs another group to transfer a member
                # Base this on probability
                ix = probs[:,k].argsort()[::-1]
                for i in ix:
                    # i is our candidate for transfer
                    # but first check that moving it
                    # doesn't empty its current group
                    curr_grp = new_assignment[i]
                    curr_grp_count = np.count_nonzero(new_assignment==curr_grp)
                    if curr_grp_count < 2:
                        logger.info('Transferring item {} would cause group {} to become empty!'.format(i, curr_grp))
                    if curr_grp_count > 1:
                        new_assignment[i] = k
                        logger.info('Transferred item {} to group {}'.format(i, k))
                        break
        return new_assignment

    def _fill_empty_groups_old(self, probs, assignment):
        """ Does the simple thing - if any group is empty, but needs to have at
        least one member, assign the data point with highest probability of
        membership """
        new_assignment = np.array(assignment.tolist())
        for k in range(self.numgrp):
            if np.count_nonzero(assignment==k) == 0:
                logger.info('Group {} became empty'.format(k))
                best = np.where(probs[:,k]==probs[:,k].max())[0][0]
                new_assignment[best] = k
                new_assignment = self._fill_empty_groups(probs, new_assignment)
        return new_assignment

    def wipe_partition(self, partition):
        """ Deletes analysis result of partition, e.g. so a repeat
        optimisation of the same partition can be done with a
        different model """
        for grp in partition.get_membership():
            grpid = self.scorer.get_id(grp)
            cache_dir = self.scorer.cache_dir
            prog = self.scorer.task_interface.name
            filename = os.path.join(cache_dir, '{}.{}.json'.format(grpid, prog))
            if os.path.exists(filename):
                os.unlink(filename)

    def likelihood_distance_matrix(self):
        # Assume all parameters are already updated
        dm = np.empty((len(self.scorer.collection), len(self.scorer.collection)))
        for i, gamma in enumerate(self.insts):
            for j, rec in enumerate(self.scorer.collection):
                gamma.set_tree(rec.tree)
                dm[i, j] = gamma.get_likelihood()
        scaled=(dm - np.diag(dm)[:,np.newaxis])
        return treeCl.DistanceMatrix.from_array(-0.5*(scaled+scaled.T))
