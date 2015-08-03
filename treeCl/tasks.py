from pllpy import pll
import tree_collection
import os
from abc import ABCMeta, abstractmethod, abstractproperty

from treeCl import treedist
from treeCl.tree import Tree
from treeCl.alignment import Alignment
from treeCl.utils.pll_helpers import pll_to_dict
from treeCl.parameters import Parameters
from treeCl.utils import fileIO, smooth_freqs
from treeCl.constants import PLL_RANDOM_SEED
from treeCl.wrappers.phylogenetics import FastTree, parse_fasttree_output
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

def fasttree_task(alignment_file, dna=False):
    fl = os.path.abspath(alignment_file)
    with fileIO.TempDir() as tmpd, fileIO.TempFile(tmpd) as treefile:
        fst = FastTree(verbose=False)
        cmd = '-gtr -gamma -pseudo -out {} {} {}'.format(treefile, '-nt' if dna else '', fl)
        fst(cmd, wait=True)
        with open(treefile) as treefl_handle:
            tree = treefl_handle.read().rstrip()
    result = parse_fasttree_output(fst.get_stderr())
    result['ml_tree'] = Tree(tree).as_string('newick', internal_labels=False, suppress_rooting=True).rstrip()
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


class FastTreeTaskInterface(TaskInterface):
    _name = 'FastTree'
    
    def scrape_args(self, records): 
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


class ApproxDistanceTaskInterface(TaskInterface):
    _name = 'Approx distance'

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
    _name = 'ML distance'

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
