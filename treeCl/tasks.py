from pllpy import pll
import itertools
import json
import tempfile
import tree_collection
import os
import random
from abc import ABCMeta, abstractmethod, abstractproperty

from treeCl import treedist
from treeCl.tree import Tree
from treeCl.alignment import Alignment
from treeCl.utils.pll_helpers import pll_to_dict
from treeCl.parameters import Parameters
from treeCl.utils import fileIO, smooth_freqs
from treeCl.constants import PLL_RANDOM_SEED
from treeCl.wrappers.phylogenetics import FastTree, parse_fasttree_output, Raxml, Phyml
from treeCl.parsers import RaxmlParser, PhymlParser
import logging
logger = logging.getLogger(__name__)
import numpy as np


class TaskInterface(object):
    __metaclass__ = ABCMeta
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


def eucdist_task(newick_string_a, newick_string_b, normalise):
    """
    Distributed version of tree_distance.eucdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.eucdist(tree_a, tree_b, normalise)
    except Exception as exc:
        eucdist_task.retry(exc=exc, countdown=1, max_retries=5)


def geodist_task(newick_string_a, newick_string_b, normalise):
    """
    Distributed version of tree_distance.geodist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.geodist(tree_a, tree_b, normalise)
    except Exception as exc:
        geodist_task.retry(exc=exc, countdown=1, max_retries=5)


def rfdist_task(newick_string_a, newick_string_b, normalise):
    """
    Distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.rfdist(tree_a, tree_b, normalise)
    except Exception as exc:
        rfdist_task.retry(exc=exc, countdown=1, max_retries=5)


def wrfdist_task(newick_string_a, newick_string_b, normalise):
    """
    Distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.wrfdist(tree_a, tree_b, normalise)
    except Exception as exc:
        wrfdist_task.retry(exc=exc, countdown=1, max_retries=5)


### TASKS that calculate trees
def pll_task(alignment_file, partition_string, guidetree=None, tree_search=True, threads=1, seed=PLL_RANDOM_SEED, frequencies=None,
             write_to_file=None):
    guidetree = True if guidetree is None else guidetree
    instance = pll(alignment_file, partition_string, guidetree, threads, seed)
    if frequencies is not None and len(frequencies) == instance.get_number_of_partitions():
        for i in range(len(frequencies)):
            instance.set_frequencies(frequencies[i], i, False)
    if tree_search:
        instance.optimise_tree_search(True)
    else:
        instance.optimise(True, True, True, True)
    result = pll_to_dict(instance)
    if write_to_file is not None: # attempt to write to file specified by write_to_file
        try:
            parameters = Parameters()
            parameters.construct_from_dict(result)
            with fileIO.fwriter(write_to_file, gz=True) as fl:
                parameters.write(fl)
        except:
            pass  # fail silently
    return result

def phyml_task(alignment_file, model):
    import re
    fl = os.path.abspath(alignment_file)
    ph = Phyml(verbose=False)
    if model in ['JC69', 'K80', 'F81', 'F84', 'HKY85', 'TN93', 'GTR']:
        datatype = 'nt'
    elif re.search('[01]{6}', model) is not None:
        datatype = 'nt'
    else:
        datatype = 'aa'
    cmd = '-i {} -m {} -d {} -f m'.format(alignment_file, model, datatype)
    ph(cmd, wait=True)
    parser = PhymlParser()
    expected_outfiles = ['{}_phyml_stats'.format(alignment_file), '{}_phyml_tree'.format(alignment_file)]
    with fileIO.TempFileList(expected_outfiles):
        try:
            result = parser.to_dict(*expected_outfiles)
        except IOError as ioerr:
            logger.error('File IO error', ioerr)
            result = None
        except ParseException as parseerr:
            logger.error('Other parse error', parseerr)
            result = None
    return result

def fasttree_task(alignment_file, dna=False):
    fl = os.path.abspath(alignment_file)
    fst = FastTree(verbose=False)
    cmd = '-gtr -gamma -pseudo {} {}'.format('-nt' if dna else '', fl)
    logger.debug('{} {}'.format(fst.exe, cmd))
    fst(cmd, wait=True)
    tree = fst.get_stdout()
    result = parse_fasttree_output(fst.get_stderr())
    result['ml_tree'] = Tree(tree).as_string('newick', internal_labels=False, suppress_rooting=True).rstrip()
    return result

def raxml_task(alignment_file, model, partitions_file=None, outfile=None, threads=1, parsimony=False, fast_tree=False):
    afl = os.path.abspath(alignment_file)
    pfl = os.path.abspath(partitions_file) if partitions_file else None
    if threads > 1:
        executable = 'raxmlHPC-PTHREADS-AVX'
        cmd = '-T {} '.format(threads)
    else:
        executable = 'raxmlHPC-AVX'
        cmd = ''
    with fileIO.TempDir() as tmpd, fileIO.TempFile(tmpd) as name:
        name = os.path.basename(name)
        rax = Raxml(executable, verbose=False)
        seed=random.randint(1000, 9999)
        outdir=os.path.abspath(tmpd)
        cmd += '-m {model} -n {name} -s {seqfile} -p {seed} -O -w {outdir}'.format(
            model=model, name=name, seqfile=afl, seed=seed,
            outdir=outdir)
        if pfl:
            cmd += ' -q {}'.format(pfl)
        if fast_tree:
            parsimony = False # fast_tree takes precedence over parsimony
            cmd += ' -f E'
        elif parsimony:
            cmd += ' -y'
        rax(cmd, wait=True)
        if fast_tree:
            # Need to follow up
            cmd = '-m {model} -n modopt -s {seqfile} -p {seed} -O -w {outdir} -f e -t {outdir}/RAxML_fastTree.{name}'.format(
                model=model, seqfile=afl, seed=seed, outdir=outdir, name=name)
            logger.debug('Follow-up cmd (fast_tree) = {}'.format(cmd))
            if pfl:
                cmd += ' -q {}'.format(pfl)
            rax(cmd, wait=True)
            parser = RaxmlParser()
            result = parser.to_dict(os.path.join(tmpd, 'RAxML_info.modopt'),
                                    os.path.join(tmpd, 'RAxML_result.modopt'), dash_f_e=True)
        elif parsimony:
            # Need to follow up
            cmd = '-m {model} -n modopt -s {seqfile} -p {seed} -O -w {outdir} -f e -t {outdir}/RAxML_parsimonyTree.{name}'.format(
                model=model, seqfile=afl, seed=seed, outdir=outdir, name=name)
            logger.debug('Follow-up cmd (parsimony) = {}'.format(cmd))
            if pfl:
                cmd += ' -q {}'.format(pfl)
            rax(cmd, wait=True)
            parser = RaxmlParser()
            result = parser.to_dict(os.path.join(tmpd, 'RAxML_info.modopt'),
                                    os.path.join(tmpd, 'RAxML_result.modopt'), dash_f_e=True)
        else:
            logger.debug('RaxML command - {}'.format(cmd))
            parser = RaxmlParser()
            logger.debug('Temp dir exists - {} ({})'.format(os.path.exists(tmpd), os.path.abspath(tmpd)))
            logger.debug('Info file exists - {}'.format('True' if os.path.exists(os.path.join(tmpd, 'RAxML_info.{}'.format(name))) else 'False'))
            logger.debug('Tree file exists - {}'.format('True' if os.path.exists(os.path.join(tmpd, 'RAxML_bestTree.{}'.format(name))) else 'False'))
            result = parser.to_dict(os.path.join(tmpd, 'RAxML_info.{}'.format(name)),
                                    os.path.join(tmpd, 'RAxML_bestTree.{}'.format(name)))
        if outfile is not None:
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

def calc_distances_task(pll_dict, alignment_file, model=None):
    rec = Alignment(alignment_file, 'phylip', True)
    freqs = smooth_freqs(pll_dict['partitions'][0]['frequencies'])
    alpha = pll_dict['partitions'][0]['alpha']
    if model is None:
        rec.set_substitution_model('GTR' if rec.is_dna() else 'LG08+F')
    else:
        rec.set_substitution_model(model)
    rec.set_gamma_rate_model(4, alpha)
    rec.set_frequencies(freqs)
    if rec.is_dna():
        rec.set_rates(pll_dict['partitions'][0]['rates'], 'ACGT')
    rec.compute_distances()
    result = dict(distances=rec.get_distances().tolist(),
                  variances=rec.get_variances().tolist())
    pll_dict['partitions'][0].update(result)
    pll_dict['nj_tree'] = rec.get_bionj_tree()
    return pll_dict


def simulate_task(n, model, frequencies, alpha, tree, rates=None):
    rec = Alignment()
    rec.set_substitution_model(model)
    rec.set_frequencies(frequencies)
    rec.set_gamma_rate_model(4, alpha)
    if rates is not None:
        try:
            rec.set_rates(rates, 'acgt')
        except RuntimeError:
            pass
    rec.set_simulator(tree)
    return rec.simulate(n)


def minsq_task(dv, gm, lab, tree, niters=10, keep_topology=False):
    tree, lk = tree_collection.compute(dv, gm, lab, tree, niters, True, keep_topology, False)
    tree = Tree(tree)
    tree.deroot()
    return dict(tree=tree.newick, score=lk)


class PllTaskInterface(TaskInterface):
    _name = 'PLL'

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
                curr_args = (filename, partition, tree, tree_search, 1, PLL_RANDOM_SEED, None, output_file)
            else:
                curr_args = (filename, partition, tree, tree_search, 1, PLL_RANDOM_SEED)
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


class RaxmlTaskInterface(TaskInterface):
    _name = 'Raxml'

    def scrape_args(self, records, partition_files=None, model=None, outfiles=None, threads=1, parsimony=False, fast_tree=False):
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
                    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
                        qfile = tmpfile.name
                        to_delete.append(tmpfile.name)
                        mymodel = 'DNAX' if rec.is_dna() else model.replace('PROT', '').replace('GAMMA', '')
                        partition_string = '{model}, {name} = 1-{seqlen}\n'.format(
                            model=mymodel,
                            name=rec.name, seqlen=len(rec))
                        tmpfile.write(partition_string)

            args.append((filename, model, qfile, ofile, threads, parsimony, fast_tree))
        return args, to_delete

    def get_task(self):
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
    """
    BROKEN
    """
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

class TreeDistanceTaskInterface(TaskInterface):
    __metaclass__ = ABCMeta
    _name = 'TreeDistance'

    def scrape_args(self, trees, normalise=False):
        return [(t1, t2, normalise) for (t1, t2) in itertools.combinations(trees, 2)]

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

