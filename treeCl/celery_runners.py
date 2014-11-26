""" definitions put here to avoid circular import between treedist and treeCl_tasks
"""
import functools
import itertools
import math
import scipy.spatial
import time
from celery import group
from treeCl.tasks.tasks import eucdist_task, geodist_task, rfdist_task, wrfdist_task
from treeCl.utils import flatten_list, setup_progressbar

def _generic_matrix_task(task, trees, normalise):

    jobs = itertools.combinations(trees, 2)
    n_jobs = int(math.ceil(len(trees) * (len(trees)-1) / 200))
    pbar = setup_progressbar("Getting inter-tree distances", n_jobs)
    pbar.start()
    #job_group = group(task.s(t1, t2, normalise) for (t1, t2) in jobs)()

    # Split the work into batches of 100. Each batch is executed asynchronously
    # -- this should work better for large amounts of quick tasks
    job_chunks = task.chunks(((t1, t2, normalise) for (t1, t2) in jobs), 100).group()()
    while not job_chunks.ready():
        pbar.update(job_chunks.completed_count())
        time.sleep(2)
    pbar.finish()
    results = job_chunks.get()
    return scipy.spatial.distance.squareform(flatten_list(results))

eucdist_matrix_task = functools.partial(_generic_matrix_task, eucdist_task)
geodist_matrix_task = functools.partial(_generic_matrix_task, geodist_task)
rfdist_matrix_task  = functools.partial(_generic_matrix_task, rfdist_task)
wrfdist_matrix_task = functools.partial(_generic_matrix_task, wrfdist_task)

