#!/usr/bin/env python
from __future__ import print_function

# standard lib
import sys
import itertools
import timeit

# third party
from dendropy import TaxonSet

# treeCl
from datastructs.trcl_seq import TrClSeq
from datastructs.trcl_tree import TrClTree
from distance_matrix import DistanceMatrix
from software_interfaces.alf import lsf_simulate_from_record, \
    simulate_from_record
from software_interfaces.DVscript import runDV
from software_interfaces.phyml import runPhyml, runLSFPhyml
from software_interfaces.treecollection import runTC
from utils import flatten_list
from utils.lazyprop import lazyprop
from errors import OptionError, optioncheck, directorymake, \
    directorycheck, isnumbercheck
from utils import fileIO
from constants import TMPDIR, SORT_KEY, ANALYSES, PHYML_MEMORY_MIN


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
            datatype=None,
            tmpdir=TMPDIR,
            calc_distances=False,
            compression=None,
            debug=False,
    ):

        self.tmpdir = directorymake(tmpdir)
        self._records = None
        self.debug = debug

        if records:
            self.records = records
            self.datatype = datatype or records[0].datatype
            optioncheck(self.datatype, ['dna', 'protein'])
            for rec in records:
                rec.tmpdir = self.tmpdir

        elif input_dir:
            directorycheck(input_dir)
            self.datatype = optioncheck(datatype, ['dna', 'protein'])
            optioncheck(file_format, ['fasta', 'phylip'])
            self.records = self.read_alignments(input_dir,
                                                file_format,
                                                compression)

        else:
            raise Exception('Provide a list of records, '
                            'or the path to a set of alignments')

        self.taxon_set = TaxonSet()
        if trees_dir:
            self.read_trees(trees_dir, self.taxon_set)

        if not self.records:
            raise NoRecordsError(file_format, input_dir, compression)

        if calc_distances:
            self.calc_distances()


    def __len__(self):
        if hasattr(self, 'records'):
            return len(self.records)
        return 0

    def __getitem__(self, i):
        if hasattr(self, 'records'):
            return self.records[i]

    @property
    def debug(self):
        return (self._debug if hasattr(self, '_debug') else False)

    @debug.setter
    def debug(self, boolean):
        self._debug = bool(boolean)

    @property
    def records(self):
        """ Returns a list of records in SORT_KEY order """
        return [self._records[i] for i in range(len(self._records))]

    @records.setter
    def records(self, records):
        """ Sets a dictionary of records keyed by SORT_KEY order """
        for rec in records:
            rec.sanitise()
        self._records = dict(enumerate(records))

        # self.reverse_lookup = {v:k for (k,v) in self.records}

    @property
    def trees(self):
        """ Returns a list of trees in SORT_KEY order """
        return [self._records[i].tree for i in range(len(self._records))]

    def num_species(self):
        """ Returns the number of species found over all records
        """
        all_headers = reduce(lambda x, y: set(x) | set(y),
                             (rec.headers for rec in self.records))
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

        return [TrClSeq(f, file_format=file_format, datatype=self.datatype,
                        name=fileIO.strip_extensions(f),
                        tmpdir=self.tmpdir)
                for f in files]

    def read_trees(self, input_dir, taxon_set=None):
        """ Read a directory full of tree files, matching them up to the
        already loaded alignments """

        extensions = ['nwk', 'tree']
        files = fileIO.glob_by_extensions(input_dir, extensions)
        files.sort(key=SORT_KEY)

        trees = [TrClTree.read_from_file(file_, taxon_set) for file_ in files]

        for record, tree in zip(self.records, trees):
            assert record.name == tree.name
            record.tree = tree

    def calc_distances(self, verbosity=0):
        """ Calculates within-alignment pairwise distances for every
        alignment. Uses Darwin. """
        for rec in self.records:
            runDV(rec, verbosity=verbosity)

    def calc_TC_trees(self, verbosity=0):
        """ Calculates distances trees using TreeCollection """
        self.analysis = 'TreeCollection'
        for rec in self.records:
            runTC(rec, self.tmpdir, verbosity=verbosity,
                  taxon_set=self.taxon_set)
            rec.tree = TrClTree.cast(rec.tree)

    def calc_TC_trees_new(self, verbosity=0):
        for rec in self.records:
            rec.tree_collection(taxon_set=self.taxon_set,
                                quiet=(True if verbosity == 0 else False))

    def calc_phyml_trees(self, analysis='nj', lsf=False, strategy='dynamic',
                         minmem=256, bootstraps=None, add_originals=False,
                         verbosity=0):
        """ Calculates trees for each record using phyml """
        optioncheck(analysis, ANALYSES)
        if bootstraps is not None:
            bootstraps = int(isnumbercheck(bootstraps))
            records = list(itertools.chain(*[[r.bootstrap_sample(str(i))
                                              for i in range(bootstraps)]
                                             for r in self]))
            if add_originals:
                records.extend(self.records)
        else:
            records = self.records

        if lsf:
            trees = runLSFPhyml(records,
                                self.tmpdir,
                                analysis=analysis,
                                verbosity=verbosity,
                                strategy=strategy,
                                minmem=minmem,
                                taxon_set=self.taxon_set,
                                debug=self.debug)
            for rec, tree in zip(records, trees):
                rec.tree = TrClTree.cast(tree)

        else:
            for rec in records:
                runPhyml(rec, self.tmpdir, analysis=analysis,
                         verbosity=verbosity, taxon_set=self.taxon_set)
                rec.tree = TrClTree.cast(rec.tree)
        if verbosity == 1:
            print()

        if bootstraps is not None:
            return [r.tree for r in records]

    def get_phyml_command_strings(self, analysis, tmpdir, verbosity=0):
        """ Gets command lines required for running phyml on every record """
        cmds = [runPhyml(rec, tmpdir, analysis=analysis,
                         verbosity=verbosity, taxon_set=self.taxon_set,
                         dry_run=True)
                for rec in self.records]
        return cmds

    def distance_matrix(self, metric, **kwargs):
        """ Generate a distance matrix from a fully-populated Collection """
        return DistanceMatrix(self.trees, metric, tmpdir=self.tmpdir,
                              **kwargs)

    def permuted_copy(self):
        """ Return a copy of the collection with all alignment columns permuted
        """
        lengths, names = zip(*[(rec.seqlength, rec.name)
                               for rec in self.records])
        concat = Concatenation(self, range(len(self))).sequence_record
        concat.shuffle()
        new_records = concat.split_by_lengths(lengths, names)
        return self.__class__(new_records)


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
        return list(itertools.chain(*[self.collection.records[i].dv
                                      for i in self.indices]))

    @lazyprop
    def sequence_record(self):
        seq0 = self.collection.records[self.indices[0]]
        for i in self.indices[1:]:
            seq0 += self.collection.records[i]
        seq0.name = '-'.join(str(x) for x in self.indices)
        return seq0

    @lazyprop
    def names(self):
        return [self.collection.records[i].name for i in self.indices]

    @lazyprop
    def lengths(self):
        return [self.collection.records[i].seqlength for i in self.indices]

    @lazyprop
    def headers(self):
        return [self.collection.records[i].headers for i in self.indices]

    @lazyprop
    def coverage(self):
        total = float(self.collection.num_species())
        return [len(self.collection.records[i]) / total for i in self.indices]

    @lazyprop
    def datatypes(self):
        return [self.collection.records[i].datatype for i in self.indices]

    def qfile(self, dna_model='GTRGAMMA', protein_model='PROTGAMMAWAG',
              ml_freqs=False):
        from_ = 1
        to_ = 0
        qs = list()
        if ml_freqs:
            dna_model += 'X'
            protein_model += 'X'

        models = dict(dna=dna_model, protein=protein_model)
        for length, name, datatype in zip(self.lengths, self.names,
                                          self.datatypes):
            to_ += length
            qs.append('{}, {} = {}-{}'.format(models[datatype], name, from_,
                                              to_))
            from_ += length
        return '\n'.join(qs)

    def paml_partitions(self):
        return 'G {} {}'.format(len(self.lengths),
                                ' '.join(str(x) for x in self.lengths))


class Scorer(object):
    """ Takes an index list, generates a concatenated SequenceRecord, calculates
    a tree and score """

    def __init__(
            self,
            collection,
            analysis,
            lsf=False,
            max_guidetrees=10,
            tmpdir=None,
            datatype=None,
            verbosity=0,
            populate_cache=True,
            debug=False,
    ):

        optioncheck(analysis, ANALYSES + ['tc', 'TreeCollection'])
        if analysis == 'tc':
            self.analysis = 'TreeCollection'
        else:
            self.analysis = analysis
        self.max_guidetrees = max_guidetrees
        self.lsf = lsf
        self.collection = collection
        self.datatype = datatype or collection.datatype
        self.verbosity = verbosity
        optioncheck(self.datatype, ['protein', 'dna'])
        self.tmpdir = tmpdir or collection.tmpdir
        directorymake(self.tmpdir)
        self.cache = {}
        self.history = []
        self.debug = debug
        if populate_cache:
            self.populate_cache()

    @property
    def records(self):
        return self.collection.records

    @property
    def debug(self):
        return (self._debug if hasattr(self, '_debug') else False)

    @debug.setter
    def debug(self, boolean):
        self._debug = bool(boolean)

    def add_partition_list(self, partition_list):
        """ Calculates concatenated trees for a list of Partitions """
        index_tuples = list(itertools.chain(*[partition.get_membership()
                                              for partition in partition_list]))
        missing = sorted(set(index_tuples).difference(self.cache.keys()))
        self._add_index_tuple_list(missing)

    def _add_index_tuple_list(self, index_tuple_list):
        if self.lsf and not self.analysis == 'TreeCollection':
            supermatrices = [self.concatenate(index_tuple).sequence_record
                             for index_tuple in index_tuple_list]

            trees = [TrClTree.cast(tree) for tree in runLSFPhyml(supermatrices,
                                                                 self.tmpdir,
                                                                 analysis=self.analysis,
                                                                 verbosity=self.verbosity,
                                                                 strategy='dynamic',
                                                                 minmem=PHYML_MEMORY_MIN,
                                                                 debug=self.debug,
                                                                 taxon_set=self.collection.taxon_set)]
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
            sequence_record = self.collection[index_tuple[0]]
        else:
            concat = self.concatenate(index_tuple)
            sequence_record = concat.sequence_record

        if self.analysis == 'TreeCollection':
            guidetrees = [self.records[n].tree for n in
                          index_tuple][:self.max_guidetrees]
            tree = sequence_record.tree_collection(
                taxon_set=self.collection.taxon_set,
                guide_tree=guidetrees[0])

        else:
            tree = TrClTree.cast(runPhyml(sequence_record,
                                          self.tmpdir,
                                          analysis=self.analysis,
                                          verbosity=self.verbosity,
                                          taxon_set=self.collection.taxon_set))
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
        concat = self.concatenate(index_tuple).sequence_record
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

    def dump(self, filename):
        """ Convenience wrapper to pickle the object. Gzipped pickle
        is written to filename """
        fileIO.gpickle(self, filename)

    def dump_cache(self, filename):
        with open(filename, 'w') as outfile:
            for k, v in self.cache.items():
                outfile.write('{}\t{}\n'.format(k, str(v)))

    def load_cache(self, filename):
        d = {}
        with open(filename, 'r') as infile:
            for line in infile:
                k, s = line.rstrip().split('\t')
                t = TrClTree.gen_from_text(s, self.collection.taxon_set)
                d[k] = t
        self.cache = d

