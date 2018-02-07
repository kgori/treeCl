from __future__ import absolute_import
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
import itertools
import json
import tempfile
import tree_collection
import os
import random
from abc import ABCMeta, abstractmethod, abstractproperty
from functools import reduce

import phylo_utils

from . import treedist
from .tree import Tree
from .alignment import Alignment, SequenceSimulator
from .parameters import Parameters
from .utils import fileIO, smooth_freqs
from .constants import RANDOM_SEED
from .wrappers.phylogenetics import FastTree, parse_fasttree_output, Raxml, Phyml
from .parsers import RaxmlParser, PhymlParser
import logging
from future.utils import with_metaclass
logger = logging.getLogger(__name__)
import numpy as np


class TaskInterface(with_metaclass(ABCMeta, object)):
    _name = None

    @abstractmethod
    def get_task(self):
        pass

    @abstractmethod
    def scrape_args(self):
        pass

    @property
    def name(self):
        return self._name

def eucdist_task(newick_string_a, newick_string_b, normalise, min_overlap=4):
    """
    Distributed version of tree_distance.eucdist
    Parameters: two valid newick strings and a boolean
    """
    tree_a = Tree(newick_string_a)
    tree_b = Tree(newick_string_b)
    return treedist.eucdist(tree_a, tree_b, normalise, min_overlap)

def geodist_task(newick_string_a, newick_string_b, normalise, min_overlap=4):
    """
    Distributed version of tree_distance.geodist
    Parameters: two valid newick strings and a boolean
    """
    tree_a = Tree(newick_string_a)
    tree_b = Tree(newick_string_b)
    return treedist.geodist(tree_a, tree_b, normalise, min_overlap)

def rfdist_task(newick_string_a, newick_string_b, normalise, min_overlap=4):
    """
    Distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    tree_a = Tree(newick_string_a)
    tree_b = Tree(newick_string_b)
    return treedist.rfdist(tree_a, tree_b, normalise, min_overlap)

def wrfdist_task(newick_string_a, newick_string_b, normalise, min_overlap=4):
    """
    Distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    tree_a = Tree(newick_string_a)
    tree_b = Tree(newick_string_b)
    return treedist.wrfdist(tree_a, tree_b, normalise, min_overlap)

### TASKS that calculate trees
def pll_task(alignment_file, partition_string, guidetree=None, tree_search=True, threads=1, seed=RANDOM_SEED, frequencies=None,
             write_to_file=None):
    try:
        import pllpy
    except:
        logger.error("Couldn't import pllpy: returning empty dict")
        retval = {}
        return retval
    guidetree = True if guidetree is None else guidetree
    instance = pllpy.pll(alignment_file, partition_string, guidetree, threads, seed)
    if frequencies is not None and len(frequencies) == instance.get_number_of_partitions():
        for i in range(len(frequencies)):
            instance.set_frequencies(frequencies[i], i, False)
    if tree_search:
        instance.optimise_tree_search(True)
    else:
        instance.optimise(True, True, True, True)
    result = pllpy.helpers.pll_to_dict(instance)
    if write_to_file is not None: # attempt to write to file specified by write_to_file
        try:
            parameters = Parameters()
            parameters.construct_from_dict(result)
            with fileIO.fwriter(write_to_file, gz=True) as fl:
                parameters.write(fl)
        except:
            pass  # fail silently
    return result

def phyml_task(alignment_file, model, **kwargs):
    """
    Kwargs are passed to the Phyml process command line
    """
    import re
    fl = os.path.abspath(alignment_file)
    ph = Phyml(verbose=False)
    if model in ['JC69', 'K80', 'F81', 'F84', 'HKY85', 'TN93', 'GTR']:
        datatype = 'nt'
    elif re.search('[01]{6}', model) is not None:
        datatype = 'nt'
    else:
        datatype = 'aa'
    cmd = '-i {} -m {} -d {} -f m --quiet'.format(alignment_file, model, datatype)
    logger.debug("Phyml command = {}".format(cmd))
    ph(cmd, wait=True, **kwargs)
    logger.debug("Phyml stdout = {}".format(ph.get_stdout()))
    logger.debug("Phyml stderr = {}".format(ph.get_stderr()))
    parser = PhymlParser()
    expected_outfiles = ['{}_phyml_stats'.format(alignment_file), '{}_phyml_tree'.format(alignment_file)]
    for i in range(2):
        if not os.path.exists(expected_outfiles[i]):
            expected_outfiles[i] += '.txt'
    logger.debug('Stats file {} {}'.format(expected_outfiles[0], 'exists' if os.path.exists(expected_outfiles[0]) else 'doesn\'t exist'))
    logger.debug('Tree file {} {}'.format(expected_outfiles[1], 'exists' if os.path.exists(expected_outfiles[1]) else 'doesn\'t exist'))
    with fileIO.TempFileList(expected_outfiles):
        try:
            result = parser.to_dict(*expected_outfiles)
        except IOError as ioerr:
            logger.error('File IO error: {}'.format(ioerr))
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error: {}'.format(parseerr))
            result = None
    return result

def bionj_task(alignment_file, model, **kwargs):
    """
    Kwargs are passed to the Phyml process command line
    """
    import re
    fl = os.path.abspath(alignment_file)
    ph = Phyml(verbose=False)
    if model in ['JC69', 'K80', 'F81', 'F84', 'HKY85', 'TN93', 'GTR']:
        datatype = 'nt'
    elif re.search('[01]{6}', model) is not None:
        datatype = 'nt'
    else:
        datatype = 'aa'
    cmd = '-i {} -m {} -d {} -b 0 -o n --quiet'.format(alignment_file, model, datatype)
    logger.debug("Phyml command = {}".format(cmd))
    ph(cmd, wait=True, **kwargs)
    logger.debug("Phyml stdout = {}".format(ph.get_stdout()))
    logger.debug("Phyml stderr = {}".format(ph.get_stderr()))
    parser = PhymlParser()
    expected_outfiles = ['{}_phyml_stats'.format(alignment_file), '{}_phyml_tree'.format(alignment_file)]
    for i in range(2):
        if not os.path.exists(expected_outfiles[i]):
            expected_outfiles[i] += '.txt'
    logger.debug('Stats file {} {}'.format(expected_outfiles[0], 'exists' if os.path.exists(expected_outfiles[0]) else 'doesn\'t exist'))
    logger.debug('Tree file {} {}'.format(expected_outfiles[1], 'exists' if os.path.exists(expected_outfiles[1]) else 'doesn\'t exist'))
    with fileIO.TempFileList(expected_outfiles):
        try:
            result = parser.to_dict(*expected_outfiles)
        except IOError as ioerr:
            logger.error('File IO error: {}'.format(ioerr))
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error: {}'.format(parseerr))
            result = None
    return result

def fasttree_task(alignment_file, dna=False):
    fl = os.path.abspath(alignment_file)
    fst = FastTree(verbose=False)
    cmd = '{} -gamma -pseudo {} {}'.format('-gtr' if dna else '-wag', '-nt' if dna else '', fl)
    logger.debug('{} {}'.format(fst.exe, cmd))
    fst(cmd, wait=True)
    tree = fst.get_stdout()
    result = parse_fasttree_output(fst.get_stderr())
    result['ml_tree'] = Tree(tree).newick
    result['partitions'][0]['model'] = 'GTR' if dna else 'WAG'
    return result

def raxml_task(executable, alignment_file, model, partitions_file=None, outfile=None, threads=1, parsimony=False, fast_tree=False):
    logger.debug('raxml_task: executable {}, alignment_file {}, model {}, partitions_file {}, outfile {}, threads {}, parsimony {}, fast_tree {}'.format(executable, alignment_file, model, partitions_file, outfile, threads, parsimony, fast_tree))
    afl = os.path.abspath(alignment_file)
    pfl = os.path.abspath(partitions_file) if partitions_file else None
    if threads > 1:
        if 'raxmlHPC' in executable and not 'PTHREADS' in executable:
            executable = executable.replace('raxmlHPC', 'raxmlHPC-PTHREADS')
            logger.debug("Sequential executable modified because threading requested: {}".format(executable))
        basecmd = '-T {} '.format(threads)
    else:
        basecmd = ''

    # initialise RAxML wrapper
    rax = Raxml(executable, verbose=False)

    with fileIO.TempDir(disable_delete=False) as tmpd, fileIO.TempFile(tmpd, disable_delete=False) as name:
        name = os.path.basename(name)
        seed=random.randint(1000, 9999)
        outdir=os.path.abspath(tmpd)
        if not os.path.exists(outdir):
            outdir = os.path.abspath('.')
        logger.debug('Raxml output files will be written to {}'.format(outdir))
        cmd = basecmd + '-m {model} -n {name} -s {seqfile} -p {seed} -O -w {outdir}'.format(
            model=model, name=name, seqfile=afl, seed=seed,
            outdir=outdir)
        if pfl:
            cmd += ' -q {}'.format(pfl)
        if fast_tree:
            parsimony = False # fast_tree takes precedence over parsimony
            cmd += ' -f E'
        elif parsimony:
            cmd += ' -y'
        logger.debug('Launching {} {}'.format(executable, cmd))
        rax(cmd, wait=True)
        logger.debug(rax.get_stdout())
        logger.debug(rax.get_stderr())

        if fast_tree:
            # Need to follow up
            fast_tree_file = os.path.join(outdir, 'RAxML_fastTree.{}'.format(name))
            logger.debug('Fast tree file exists {}'.format('yes' if os.path.exists(fast_tree_file) else 'no'))
            cmd = basecmd + '-m {model} -n modopt -s {seqfile} -p {seed} -O -w {outdir} -f e -t {fast_tree_file}'.format(
                model=model, seqfile=afl, seed=seed, outdir=outdir, fast_tree_file=fast_tree_file)
            if pfl:
                cmd += ' -q {}'.format(pfl)
            logger.debug('Launching fast tree follow-up: {} {}'.format(executable, cmd))
            rax(cmd, wait=True)
            logger.debug(rax.get_stdout())
            logger.debug(rax.get_stderr())
            info_file = os.path.join(outdir, 'RAxML_info.modopt')
            result_file = os.path.join(outdir, 'RAxML_result.modopt')
            dash_f_e = True

        elif parsimony:
            # Need to follow up
            parsimony_tree_file = os.path.join(outdir, 'RAxML_parsimonyTree.{}'.format(name))
            logger.debug('Parsimony tree file exists {}'
                         .format('yes' if os.path.exists(parsimony_tree_file) else 'no'))
            cmd = basecmd + '-m {model} -n modopt -s {seqfile} -p {seed} -O -w {outdir} -f e -t {parsimony}'.format(
                model=model, seqfile=afl, seed=seed, outdir=outdir, parsimony=parsimony_tree_file)
            if pfl:
                cmd += ' -q {}'.format(pfl)
            logger.debug('Launching parsimony follow-up: {} {}'.format(executable, cmd))
            rax(cmd, wait=True)
            logger.debug(rax.get_stdout())
            logger.debug(rax.get_stderr())
            info_file = os.path.join(outdir, 'RAxML_info.modopt')
            result_file = os.path.join(outdir, 'RAxML_result.modopt')
            dash_f_e = True

        else:
            info_file = os.path.join(outdir, 'RAxML_info.{}'.format(name))
            result_file = os.path.join(outdir, 'RAxML_result.{}'.format(name))
            dash_f_e = False

        logger.debug('Info file found - {}'.format('yes' if os.path.exists(info_file) else 'no'))
        logger.debug('Result file found - {}'.format('yes' if os.path.exists(result_file) else 'no'))
        logger.debug('Output directory found - {}'.format('yes' if os.path.isdir(outdir) else 'no'))

        if not os.path.exists(info_file):
            info_file = os.path.join(os.path.abspath('.'), 'RAxML_info.{}'.format(name))
            logger.debug('Fallback info file ({}) found - {}'
                         .format(info_file, 'yes' if os.path.exists(info_file) else 'no'))

        if not os.path.exists(result_file):
            result_file = os.path.join(os.path.abspath('.'), 'RAxML_result.{}'.format(name))
            logger.debug('Fallback result file ({}) found - {}'
                         .format(result_file, 'yes' if os.path.exists(result_file) else 'no'))

        parser = RaxmlParser()
        result = parser.to_dict(info_file, result_file, dash_f_e=dash_f_e)

        if outfile is not None:
            logger.debug('Attempting to write result to {}'.format(outfile))
            try:
                with open(outfile, 'w') as ofl:
                    json.dump(result, ofl)
            except:
                logger.error('Could not write outfile {}'.format(outfile))

    return result

def fast_calc_distances_task(alignment_file):
    rec = Alignment(alignment_file, 'phylip', True)
    rec.fast_compute_distances()
    inner = dict(distances=rec.get_distances().tolist(),
                 variances=rec.get_variances().tolist())
    outer = dict(tree=rec.get_bionj_tree(),
                 partitions={0: inner})
    return result

def calc_distances_task(parameter_dict, alignment_file, model=None):
    rec = Alignment(alignment_file, 'phylip', True)
    freqs = smooth_freqs(parameter_dict['partitions'][0]['frequencies'])
    alpha = parameter_dict['partitions'][0]['alpha']
    if model is None:
        rec.set_substitution_model('GTR' if rec.is_dna() else 'LG08+F')
    else:
        rec.set_substitution_model(model)
    rec.set_gamma_rate_model(4, alpha)
    rec.set_frequencies(freqs)
    if rec.is_dna():
        rec.set_rates(parameter_dict['partitions'][0]['rates'], 'ACGT')
    rec.compute_distances()
    result = dict(distances=rec.get_distances().tolist(),
                  variances=rec.get_variances().tolist())
    parameter_dict['partitions'][0].update(result)
    parameter_dict['nj_tree'] = rec.get_bionj_tree()
    return parameter_dict


def simulate_task(n, model, frequencies, alpha, tree, rates=None):
    if model in ('LG', 'LG08'):
        subst_model = phylo_utils.models.LG(freqs=frequencies)
    elif model == 'GTR':
        if rates is not None:
            subst_model = phylo_utils.models.GTR(freqs=frequencies, rates=rates)
        else:
            subst_model = phylo_utils.models.GTR(freqs=frequencies)

    else:
        raise ValueError('Currently only supports LG model for proteins, GTR for nucleotides')

    tmat = phylo_utils.markov.TransitionMatrix(subst_model)
    if not isinstance(tree, Tree):
        tree = Tree(tree)
    sim = SequenceSimulator(tmat, tree, ncat=4, alpha=alpha)
    return sim.simulate(n)


def minsq_task(dv, gm, lab, tree, niters=10, keep_topology=False):
    tree, lk = tree_collection.compute(dv, gm, lab, tree, niters, True, keep_topology, False)
    tree = Tree(tree)
    tree.deroot()
    return dict(tree=tree.newick, score=lk)


class PllTaskInterface(TaskInterface):
    _name = 'PLL'

    def __init__(self):
        try:
            import pllpy
        except ImportError:
            logger.error('Could not import pllpy: PllTask not available')

    def scrape_args(self, records, model=None, output_dir=None, tree_search=True, **kwargs):
        args = []
        to_delete = []
        for rec in records:
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
                curr_args = (filename, partition, tree, tree_search, 1, RANDOM_SEED, None, output_file)
            else:
                curr_args = (filename, partition, tree, tree_search, 1, RANDOM_SEED)
            args.append(curr_args)
        return args, to_delete

    def get_task(self):
        return pll_task


class PhymlTaskInterface(TaskInterface):
    _name = 'Phyml'

    def scrape_args(self, records, model=None, **kwargs):
        DEFAULT_DNA_MODEL = 'GTR'
        DEFAULT_PROTEIN_MODEL = 'LG'
        args = []
        to_delete = []
        for rec in records:
            if model is None:
                model = DEFAULT_DNA_MODEL if rec.is_dna() else DEFAULT_PROTEIN_MODEL
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            args.append((filename, model))
        return args, to_delete

    def get_task(self):
        return phyml_task


class BionjTaskInterface(TaskInterface):
    _name = 'Bionj'

    def scrape_args(self, records, model=None, **kwargs):
        DEFAULT_DNA_MODEL = 'GTR'
        DEFAULT_PROTEIN_MODEL = 'LG'
        args = []
        to_delete = []
        for rec in records:
            if model is None:
                model = DEFAULT_DNA_MODEL if rec.is_dna() else DEFAULT_PROTEIN_MODEL
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            args.append((filename, model))
        return args, to_delete

    def get_task(self):
        return bionj_task


class RaxmlTaskInterface(TaskInterface):
    """
    Provides a high-level interface to the RAxML tree inference program.
    Using the interface is a two-stage process. First, call the `scrape_args`
    method on a list of `Alignment` objects to generate a list of command line
    arguments for RAxML, and a list of temporary files that will need cleaning
    up. Second, use the `get_task` method to access a `raxml_task` function that
    can be passed, with the argument list, to a JobHandler to run the jobs.
    """
    _name = 'Raxml'

    def scrape_args(self, records, executable='raxmlHPC-AVX', partition_files=None,
                    model=None, outfiles=None, threads=1, parsimony=False, fast_tree=False):
        """
        Examine a list of records and generate RAxML command line arguments for tree inference.
        :param records: list of `Alignment` records
        :param executable: name of the RAxML executable on the system to use. Must be in the user's path.
        :param partition_files: List of RAxML partition files used to describe any partitioning scheme
            to be used (optional)
        :param model: Choice of model to use. Defaults to GTRGAMMA for DNA, or PROTGAMMALGX for amino acid alignments.
        :param outfiles: A list of output file locations to write results (required length = 1 per alignment)
        :param threads: Number of threads for RAxML to use. This is independent of any threading used by the
            `JobHandler`, and the user should be sure that their choice is appropriate for the number of threads
            available to their system, and for the RAxML executable being used.
        :param parsimony: Use RAxML's parsimony tree search only
        :param fast_tree: Use RAxML's experimental fast tree search (-f E)
        :return: (List of command line arguments, List of created temporary files)
        """
        args = []
        to_delete = []
        if partition_files is None:
            partition_files = [None for rec in records]
        if outfiles is None:
            outfiles = [None for rec in records]
        for (rec, qfile, ofile) in zip(records, partition_files, outfiles):
            if model is None:
                model = 'GTRGAMMA' if rec.is_dna() else 'PROTGAMMALGX'
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
                to_delete.append(filename + '.reduced')

            if qfile is None:
                # Attempt to find partition file on disk, using extension 'partitions.txt'
                if filename.endswith('.phy'):
                    likely_qfile = filename.replace('phy', 'partitions.txt')
                else:
                    likely_qfile = filename + '.partitions.txt'
                if os.path.exists(likely_qfile):
                    qfile = likely_qfile
                else:
                    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmpfile:
                        qfile = tmpfile.name
                        to_delete.append(tmpfile.name)
                        mymodel = 'DNAX' if rec.is_dna() else model.replace('PROT', '').replace('GAMMA', '').replace('CAT', '')
                        partition_string = '{model}, {name} = 1-{seqlen}\n'.format(
                            model=mymodel,
                            name=rec.name, seqlen=len(rec))
                        tmpfile.write(partition_string)

            args.append((executable, filename, model, qfile, ofile, threads, parsimony, fast_tree))
        return args, to_delete

    def get_task(self):
        """
        Access the RAxML command line interface wrapper
        :return: `raxml_task` function used to dispatch jobs to system RAxML
        """
        return raxml_task


class FastTreeTaskInterface(TaskInterface):
    _name = 'FastTree'

    def scrape_args(self, records, **kwargs):
        args = []
        to_delete = []
        for rec in records:
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)

            curr_args = (filename, rec.is_dna())
            args.append(curr_args)
        return args, to_delete

    def get_task(self):
        return fasttree_task


class TreeCollectionTaskInterface(TaskInterface):
    _name = 'TreeCollection'

    def scrape_args(self, records, scale=1, guide_tree=None, niters=10, keep_topology=False):
        # local lists
        distances = []
        variances = []
        headers = []
        for rec in records:
            distances.append(rec.parameters.partitions.distances)
            variances.append(rec.parameters.partitions.variances)
            headers.append(rec.get_names())

        num_matrices = len(records)
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

        tree_string = guide_tree.scale(scale).newick.replace('\'', '')

        return distvar_string, genome_map_string, labels_string, tree_string, niters, keep_topology

    def get_task(self):
        return minsq_task


class ApproxDistanceTaskInterface(TaskInterface):
    _name = 'ApproxDistance'

    def scrape_args(self, records):
        args = []
        to_delete = []
        for rec in records:
            filename, delete = rec.get_alignment_file(as_phylip=True)
            if delete:
                to_delete.append(filename)
            args.append((filename,))
        return args, to_delete

    def get_task(self):
        return fast_calc_distances_task


class MLDistanceTaskInterface(TaskInterface):
    _name = 'MLDistance'

    def scrape_args(self, records):
        args = []
        to_delete = []
        for rec in records:
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
        return args, to_delete

    def get_task(self):
        return calc_distances_task


class SimulatorTaskInterface(TaskInterface):
    _name = 'Simulator'

    def scrape_args(self, results_list):
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

class TreeDistanceTaskInterface(with_metaclass(ABCMeta, TaskInterface)):
    _name = 'TreeDistance'

    def scrape_args(self, trees, normalise=False, min_overlap=4):
        return [(t1, t2, normalise, min_overlap) for (t1, t2) in itertools.combinations(trees, 2)]

    @abstractmethod
    def get_task(self):
        pass

class GeodesicTreeDistance(TreeDistanceTaskInterface):
    _name = 'GeodesicDistance'
    def get_task(self):
        return geodist_task

class RobinsonFouldsTreeDistance(TreeDistanceTaskInterface):
    _name = 'RFDistance'
    def get_task(self):
        return rfdist_task

class WeightedRobinsonFouldsTreeDistance(TreeDistanceTaskInterface):
    _name = 'WeightedRFDistance'
    def get_task(self):
        return wrfdist_task

class EuclideanTreeDistance(TreeDistanceTaskInterface):
    _name = 'EuclideanDistance'
    def get_task(self):
        return eucdist_task


from tree_distance import getEuclideanDistance, getGeodesicDistance, getRobinsonFouldsDistance,\
    getWeightedRobinsonFouldsDistance

def _fast_geo(tree1, tree2, normalise=False):
    return getGeodesicDistance(tree1, tree2, normalise)

def _fast_euc(tree1, tree2, normalise=False):
    return getEuclideanDistance(tree1, tree2, normalise)

def _fast_rf(tree1, tree2, normalise=False):
    return getRobinsonFouldsDistance(tree1, tree2, normalise)

def _fast_wrf(tree1, tree2, normalise=False):
    return getWeightedRobinsonFouldsDistance(tree1, tree2, normalise)


class EqualLeafSetGeodesicTreeDistance(TreeDistanceTaskInterface):
    _name = 'GeodesicDistance'
    def get_task(self):
        return _fast_geo


class EqualLeafSetEuclideanTreeDistance(TreeDistanceTaskInterface):
    _name = 'EuclideanDistance'
    def get_task(self):
        return _fast_euc


class EqualLeafSetRobinsonFouldsTreeDistance(TreeDistanceTaskInterface):
    _name = 'RobinsonFouldsDistance'
    def get_task(self):
        return _fast_rf


class EqualLeafSetWeightedRobinsonFouldsTreeDistance(TreeDistanceTaskInterface):
    _name = 'WeightedRobinsonFouldsDistance'
    def get_task(self):
        return _fast_wrf
