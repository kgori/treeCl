""" definitions put here to avoid circular import between treedist and treeCl_tasks
"""
import functools
import itertools
import scipy.spatial
from celery import group
from treeCl.tasks.tasks import eucdist_task, geodist_task, rfdist_task, wrfdist_task

def _generic_matrix_task(task, trees, normalise):

    jobs = itertools.combinations((tree.newick for tree in trees), 2)
    job_group = group(task.s(t1, t2, normalise) for (t1, t2) in jobs)()
    results = job_group.get()
    return scipy.spatial.distance.squareform(results)

eucdist_matrix_task = functools.partial(_generic_matrix_task, eucdist_task)
geodist_matrix_task = functools.partial(_generic_matrix_task, geodist_task)
rfdist_matrix_task  = functools.partial(_generic_matrix_task, rfdist_task)
wrfdist_matrix_task = functools.partial(_generic_matrix_task, wrfdist_task)
