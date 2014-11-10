from pllpy import pll

from treeCl import treedist
from treeCl.tree import Tree
from treeCl.tasks.celery import app
from treeCl.interfacing.pll import pll_to_dict
from treeCl.constants import PLL_RANDOM_SEED


@app.task(time_limit=5)
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


@app.task(time_limit=5)
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


@app.task(time_limit=5)
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


@app.task(time_limit=5)
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


@app.task()
def pll_unpartitioned_task(alignment_file, partition_string, guidetree=None, threads=1, seed=PLL_RANDOM_SEED):
    guidetree = True if guidetree is None else guidetree
    instance = pll(alignment_file, partition_string, guidetree, threads, seed)
    # instance.set_optimisable_frequencies(0, True)
    instance.optimise_tree_search(True)
    return pll_to_dict(instance)
