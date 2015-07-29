from abc import ABCMeta, abstractmethod, abstractproperty
from .constants import PARALLEL_PROFILE
from .utils import setup_progressbar, grouper, flatten_list
import logging
logger = logging.getLogger(__name__)

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


def tupleise(args):
    return [a if isinstance(a, (tuple, list)) else (a,) for a in args]


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


def threadpool_map(task, args, message, concurrency, batchsize=1):
    """
    Helper to map a function over a range of inputs, using a threadpool, with a progress meter
    """
    import concurrent.futures

    njobs = len(args)
    batches = grouper(batchsize, tupleise(args))
    batched_task = lambda batch: [task(*job) for job in batch]
    PROGRESS = message is not None
    if PROGRESS:
        message += ' ({} workers)'.format(concurrency)
        pbar = setup_progressbar(message, njobs, simple_progress=True)
        pbar.start()
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = []
        for batch in batches:
            futures.append(executor.submit(batched_task, batch))

        if PROGRESS:
            for i, fut in enumerate(concurrent.futures.as_completed(futures), start=1):
                pbar.update(i*batchsize)

        else:
            concurrent.futures.wait(futures)

    if PROGRESS:
        pbar.finish()

    return flatten_list([fut.result() for fut in futures])


class JobHandler(object):
    """
    Base class to provide uniform interface for all job handlers
    """
    metaclass = ABCMeta

    # default values
    concurrency = 1
    batchsize = 1
    client = None

    def __init__(self, concurrency=None, client=None):
        if concurrency is not None:
            self.concurrency = concurrency

        if client is not None:
            self.client = client

    @abstractmethod
    def __call__(self, task, params, message, batchsize):
        """ If you define a message, then progress will be written to stderr """
        pass

class SequentialJobHandler(JobHandler):
    """
    Jobs are handled using a simple map
    """
    def __init__(self, *args, **kwargs):
        super(SequentialJobHandler, self).__init__()

    def __call__(self, task, params, message, batchsize):
        return sequential_map(task, args, message)


class ThreadpoolJobHandler(JobHandler):
    """
    Jobs are handled by a threadpool using concurrent.futures
    """
    def __init__(self, concurrency, *args, **kwargs):
        super(ThreadpoolJobHandler, self).__init__(concurrency=concurrency)

    def __call__(self, task, params, message, *args, **kwargs):
        return threadpool_map(task, args, message, self.concurrency)


class IPythonJobHandler(JobHandler):
    """
    Jobs are handled using an IPython.parallel.Client
    """
    def __init__(self, *args, **kwargs):
        client=get_client()
        if not client:
            logger.warn('Could not obtain an IPython parallel Client - this IPythonJobHandler will behave as a SequentialJobHandler')
        super(IPythonJobHandler, self).__init__(client=client)

    def __call__(self, task, params, message, *args, **kwargs):
        if self.client is None:
            return sequential_map(task, params, message, *args, **kwargs)
        else:
            return parallel_map(task, params, message, *args, **kwargs)
