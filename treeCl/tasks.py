from pllpy import pll
import tree_collection
import os

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

def eucdist_task(newick_string_a, newick_string_b, normalise):
    """
    Celery-distributed version of tree_distance.eucdist
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
    Celery-distributed version of tree_distance.geodist
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
    Celery-distributed version of tree_distance.rfdist
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
    Celery-distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.wrfdist(tree_a, tree_b, normalise)
    except Exception as exc:
        wrfdist_task.retry(exc=exc, countdown=1, max_retries=5)


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
    logger.debug("alignment_file = {}".format(alignment_file))   
    fl = os.path.abspath(alignment_file)
    logger.debug("fl = {}".format(fl))
    with fileIO.TempDir() as tmpd, fileIO.TempFile(tmpd) as treefile:
        logger.debug("Current dir = {}".format(os.getcwd()))
        fst = FastTree(verbose=False)
        cmd = '-gtr -gamma -pseudo -out {} {} {}'.format(treefile, '-nt' if dna else '', fl)
        fst(cmd, wait=True)
        logger.debug("{} exists? {}".format(treefile, 'True' if os.path.exists(treefile) else 'False'))
        logger.debug(fst.get_stderr())
        with open(treefile) as treefl_handle:
            tree = treefl_handle.read().rstrip()
        logger.debug('Tree - {}'.format(tree))
    result = parse_fasttree_output(fst.get_stderr())
    result['ml_tree'] = Tree(tree).as_string('newick', internal_labels=False, suppress_rooting=True).rstrip()
    return result


def fasttree_task_notmpdir(alignment_file, dna=False):
    fl = os.path.abspath(alignment_file)
    with fileIO.TempFile() as treefile:
        fst = FastTree(verbose=False)
        cmd = '-gtr -gamma -pseudo -out {} {} {}'.format(treefile, '-nt' if dna else '', fl)
        fst(cmd, wait=True)
        with open(treefile) as treefl_handle:
            tree = treefl_handle.read().rstrip()
        logger.debug('Tree - {}'.format(tree))
    result = parse_fasttree_output(fst.get_stderr())
    result['ml_tree'] = Tree(tree).as_string('newick', internal_labels=False, suppress_rooting=True).rstrip()
    return result


def fast_calc_distances_task(alignment_file):
    rec = Alignment(alignment_file, 'phylip', True)
    rec.fast_compute_distances()
    result = dict(distances=rec.get_distances().tolist(),
                  variances=rec.get_variances().tolist(),
                  tree=rec.get_bionj_tree())
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
    tree, lk = tree_collection.compute(dv, gm, lab, tree, niters, True, keep_topology, True)
    tree = Tree(tree)
    tree.deroot()
    return dict(tree=tree.newick, score=lk)
