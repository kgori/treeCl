""" definitions put here to avoid circular import between treedist and treeCl_tasks
(This is such spaghetti code - todo: fix this shit)
"""
import functools
import itertools
import math
import scipy.spatial
import time
from treeCl.parallel.tasks import eucdist_task, geodist_task, rfdist_task, wrfdist_task
from treeCl.utils import flatten_list, setup_progressbar

__all__ = [
    'eucdist_matrix_async',
    'geodist_matrix_async',
    'rfdist_matrix_async',
    'wrfdist_matrix_async',
    'eucdist_matrix_sequential',
    'geodist_matrix_sequential',
    'rfdist_matrix_sequential',
    'wrfdist_matrix_sequential',
]

def _generic_async_matrix_task(task, trees, normalise, batch_size=100):

    jobs = itertools.combinations(trees, 2)
    n_jobs = int(math.ceil(len(trees) * (len(trees)-1) / (2*batch_size)))
    pbar = setup_progressbar("Getting inter-tree distances (async)", n_jobs, simple_progress=True)
    pbar.start()

    # Split the work into batches of 'batch_size'. Each batch is executed asynchronously
    # -- this should work better for large amounts of quick tasks
    job_chunks = task.chunks(((t1, t2, normalise) for (t1, t2) in jobs), batch_size).group()()
    while not job_chunks.ready():
        pbar.update(job_chunks.completed_count())
        time.sleep(2)
    pbar.finish()
    results = job_chunks.get()
    return scipy.spatial.distance.squareform(flatten_list(results))

def _generic_sequential_matrix_task(task, trees, normalise):
    jobs = itertools.combinations(trees, 2)
    n_jobs = int(math.ceil(len(trees) * (len(trees)-1) / 2))
    pbar = setup_progressbar("Getting inter-tree distances (seq)", n_jobs)
    pbar.start()
    results = []
    for i, (t1, t2, normalise) in enumerate((t1, t2, normalise) for (t1, t2) in jobs):
        results.append(task(t1, t2, normalise))
        pbar.update(i)

    pbar.finish()
    return scipy.spatial.distance.squareform(results)

eucdist_matrix_async = functools.partial(_generic_async_matrix_task, eucdist_task)
geodist_matrix_async = functools.partial(_generic_async_matrix_task, geodist_task)
rfdist_matrix_async  = functools.partial(_generic_async_matrix_task, rfdist_task)
wrfdist_matrix_async = functools.partial(_generic_async_matrix_task, wrfdist_task)
eucdist_matrix_sequential = functools.partial(_generic_sequential_matrix_task, eucdist_task)
geodist_matrix_sequential = functools.partial(_generic_sequential_matrix_task, geodist_task)
rfdist_matrix_sequential  = functools.partial(_generic_sequential_matrix_task, rfdist_task)
wrfdist_matrix_sequential = functools.partial(_generic_sequential_matrix_task, wrfdist_task)
