from treeCl.constants import PARALLEL_PROFILE
from treeCl.utils import setup_progressbar

__author__ = 'kgori'


def async_avail():
    from IPython import parallel
    try:
        client = parallel.Client(PARALLEL_PROFILE)
        return len(client) > 0
    except IOError:
        return False
    except Exception:
        return False


def get_client():
    from IPython import parallel
    try:
        client = parallel.Client(profile=PARALLEL_PROFILE)
        return client if len(client) > 0 else None
    except IOError:
        return None
    except Exception:
        return None


def parallel_map(client, task, args, message, batchsize=1, background=False):
    """
    Helper to map a function over a sequence of inputs, in parallel, with progress meter.
    :param client: IPython.parallel.Client instance
    :param task: Function
    :param args: Must be a list of tuples of arguments that the task function will be mapped onto.
                 If the function takes a single argument, it still must be a 1-tuple.
    :param message: String for progress bar
    :param batchsize: Jobs are shipped in batches of this size. Higher numbers mean less network traffic,
                      but longer execution time per job.
    :return: IPython.parallel.AsyncMapResult
    """
    njobs = len(args)
    nproc = len(client)
    view = client.load_balanced_view()
    message += ' ({} proc)'.format(nproc)
    pbar = setup_progressbar(message, njobs, simple_progress=True)
    if not background:
        pbar.start()
    map_result = view.map(task, *list(zip(*args)), chunksize=batchsize)
    if background:
        return map_result, client
    while not map_result.ready():
        map_result.wait(1)
        pbar.update(min(njobs, map_result.progress * batchsize))
    pbar.finish()
    return map_result


def sequential_map(task, args, message):
    """
    Helper to map a function over a sequence of inputs, sequentially, with progress meter.
    :param client: IPython.parallel.Client instance
    :param task: Function
    :param args: Must be a list of tuples of arguments that the task function will be mapped onto.
                 If the function takes a single argument, it still must be a 1-tuple.
    :param message: String for progress bar
    :param batchsize: Jobs are shipped in batches of this size. Higher numbers mean less network traffic,
                      but longer execution time per job.
    :return: IPython.parallel.AsyncMapResult
    """
    njobs = len(args)
    pbar = setup_progressbar(message, njobs, simple_progress=True)
    pbar.start()
    map_result = []
    for (i, arglist) in enumerate(args):
        map_result.append(task(*arglist))
        pbar.update(i)
    pbar.finish()
    return map_result
