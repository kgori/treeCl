from app import app
from celery.exceptions import SoftTimeLimitExceeded
from pllpy import pll
from treeCl import treedist, Tree
from treeCl.interfacing.pll import pll_to_dict


@app.task(time_limit=5)
def eucdist(newick_string_a, newick_string_b, normalise):
    """
    Celery-distributed version of tree_distance.eucdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.eucdist(tree_a, tree_b, normalise)
    except Exception as exc:
        eucdist.retry(exc=exc, countdown=1, max_retries=5)


@app.task(time_limit=5)
def geodist(newick_string_a, newick_string_b, normalise):
    """
    Celery-distributed version of tree_distance.geodist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.geodist(tree_a, tree_b, normalise)
    except Exception as exc:
        geodist.retry(exc=exc, countdown=1, max_retries=5)


@app.task(time_limit=5)
def rfdist(newick_string_a, newick_string_b, normalise):
    """
    Celery-distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.rfdist(tree_a, tree_b, normalise)
    except Exception as exc:
        rfdist.retry(exc=exc, countdown=1, max_retries=5)


@app.task(time_limit=5)
def wrfdist(newick_string_a, newick_string_b, normalise):
    """
    Celery-distributed version of tree_distance.rfdist
    Parameters: two valid newick strings and a boolean
    """
    try:
        tree_a = Tree(newick_string_a)
        tree_b = Tree(newick_string_b)
        return treedist.wrfdist(tree_a, tree_b, normalise)
    except Exception as exc:
        wrfdist.retry(exc=exc, countdown=1, max_retries=5)


@app.task()
def pll_unpartitioned(alignment_file, model='LG', guidetree=None):
    with open(alignment_file) as fl:
        nseq, lseq = fl.readline().split()
    partition_string = '{}, all = 1-{}\n'.format(model, lseq)
    guidetree = True if guidetree is None else guidetree
    instance = pll(alignment_file, partition_string, guidetree, 1, random.randint(10000, 100000))
    instance.set_optimisable_frequencies(0, True)
    instance.optimise(True)
    return pll_to_dict(instance)
