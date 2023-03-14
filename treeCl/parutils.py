from abc import ABCMeta, abstractmethod
from .constants import PARALLEL_PROFILE
from .utils import setup_progressbar, grouper, flatten_list
import logging
import multiprocessing
import sys
logger = logging.getLogger(__name__)

__author__ = 'kgori'

"""
Introduced this workaround for a bug in multiprocessing where
errors are thrown for an EINTR interrupt.
Workaround taken from http://stackoverflow.com/a/5395277 - but
changed because can't subclass from multiprocessing.Queue (it's
a factory method)
"""
import errno

def retry_on_eintr(function, *args, **kw):
    while True:
        try:
            return function(*args, **kw)
        except IOError as e:
            if e.errno == errno.EINTR:
                continue
            else:
                raise

def get_from_queue(queue, block=True, timeout=None):
    return retry_on_eintr(queue.get, block, timeout)
"""
End of workaround
"""

def fun(f, q_in, q_out):
    while True:
        (i, x) = get_from_queue(q_in)
        if i is None:
            break
        q_out.put((i, f(*x)))

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
    for a in args:
        if isinstance(a, (tuple, list)):
            yield a
        else:
            yield (a,)

def get_njobs(nargs, args):
    if nargs is not None:
        njobs = nargs
    elif isinstance(args, (tuple, list)):
        njobs = len(args)
    else:
        njobs = int(sys.maxsize / 1000000)  # sys.maxsize is too large for progressbar to display ETA (datetime issue)
    return njobs


def parallel_map(client, task, args, message, batchsize=1, background=False, nargs=None):
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
    show_progress = bool(message)
    njobs = get_njobs(nargs, args)
    nproc = len(client)
    logger.debug('parallel_map: len(client) = {}'.format(len(client)))
    view = client.load_balanced_view()
    if show_progress:
        message += ' (IP:{}w:{}b)'.format(nproc, batchsize)
        pbar = setup_progressbar(message, njobs, simple_progress=True)
        if not background:
            pbar.start()
    map_result = view.map(task, *list(zip(*args)), chunksize=batchsize)
    if background:
        return map_result, client
    while not map_result.ready():
        map_result.wait(1)
        if show_progress:
            pbar.update(min(njobs, map_result.progress * batchsize))
    if show_progress:
        pbar.finish()
    return map_result

def sequential_map(task, args, message, nargs=None):
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
    njobs = get_njobs(nargs, args)
    show_progress = bool(message)
    if show_progress:
        pbar = setup_progressbar(message, njobs, simple_progress=True)
        pbar.start()
    map_result = []
    for (i, arglist) in enumerate(tupleise(args), start=1):
        map_result.append(task(*arglist))
        if show_progress:
            pbar.update(i)
    if show_progress:
        pbar.finish()
    return map_result

def threadpool_map(task, args, message, concurrency, batchsize=1, nargs=None):
    """
    Helper to map a function over a range of inputs, using a threadpool, with a progress meter
    """
    import concurrent.futures


    njobs = get_njobs(nargs, args)
    show_progress = bool(message)
    batches = grouper(batchsize, tupleise(args))
    batched_task = lambda batch: [task(*job) for job in batch]
    if show_progress:
        message += ' (TP:{}w:{}b)'.format(concurrency, batchsize)
        pbar = setup_progressbar(message, njobs, simple_progress=True)
        pbar.start()
    with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = []
        completed_count = 0
        for batch in batches:
            futures.append(executor.submit(batched_task, batch))

        if show_progress:
            for i, fut in enumerate(concurrent.futures.as_completed(futures), start=1):
                completed_count += len(fut.result())
                pbar.update(completed_count)

        else:
            concurrent.futures.wait(futures)

    if show_progress:
        pbar.finish()

    return flatten_list([fut.result() for fut in futures])

def processpool_map(task, args, message, concurrency, batchsize=1, nargs=None):
    """
    See http://stackoverflow.com/a/16071616
    """
    njobs = get_njobs(nargs, args)
    show_progress = bool(message)
    batches = grouper(batchsize, tupleise(args))

    # batched_task and multiprocessing ForkingPickler don't play nicely together
    # def batched_task(*batch):
    #     return [task(*job) for job in batch]

    if show_progress:
        message += ' (PP:{}w:{}b)'.format(concurrency, batchsize)
        pbar = setup_progressbar(message, njobs, simple_progress=True)
        pbar.start()
    
    q_in   = multiprocessing.Queue()  # Should I limit either queue size? Limiting in-queue
    q_out  = multiprocessing.Queue()  # increases time taken to send jobs, makes pbar less useful

    # Not using batched_task here because it doesn't play nicely with ForkingPickler
    proc = [multiprocessing.Process(target=fun, args=(task, q_in, q_out)) for _ in range(concurrency)]
    for p in proc:
        p.daemon = True
        p.start()
    sent = [q_in.put((i, x)) for (i, x) in enumerate(args)]
    [q_in.put((None, None)) for _ in range(concurrency)]
    res = []
    completed_count = 0
    for _ in range(len(sent)):
        result = get_from_queue(q_out)
        res.append(result)
        completed_count += 1 #len(result[1])
        if show_progress:
            pbar.update(completed_count)

    [p.join() for p in proc]
    if show_progress:
        pbar.finish()

    return flatten_list([x for (i, x) in sorted(res)])


class JobHandler(object):
    """
    Base class to provide uniform interface for all job handlers
    """
    metaclass = ABCMeta

    @abstractmethod
    def __call__(self, task, args, message, batchsize):
        """ If you define a message, then progress will be written to stderr """
        pass


class SequentialJobHandler(JobHandler):
    """
    Jobs are handled using a simple map
    """
    def __call__(self, task, args, message, batchsize, nargs=None):
        if batchsize > 1:
            logger.warn("Setting batchsize > 1 has no effect when using a SequentialJobHandler")
        return sequential_map(task, args, message, nargs)


class ThreadpoolJobHandler(JobHandler):
    """
    Jobs are handled by a threadpool using concurrent.futures
    """
    def __init__(self, concurrency):
        self.concurrency = concurrency

    def __call__(self, task, args, message, batchsize, nargs=None):
        return threadpool_map(task, args, message, self.concurrency, batchsize, nargs)


class ProcesspoolJobHandler(JobHandler):
    """
    Jobs are handled by a threadpool using concurrent.futures
    """
    def __init__(self, concurrency):
        self.concurrency = concurrency

    def __call__(self, task, args, message, batchsize, nargs=None):
        return processpool_map(task, args, message, self.concurrency, batchsize, nargs)


class IPythonJobHandler(JobHandler):
    """
    Jobs are handled using an IPython.parallel.Client
    """
    def __init__(self, profile=None):
        """
        Initialise the IPythonJobHandler using the given ipython profile.

        Parameters
        ----------

        profile:  string
                  The ipython profile to connect to - this should already be running an ipcluster
                  If the connection fails it raises a RuntimeError
        """
        import IPython.parallel
        try:
            self.client=IPython.parallel.Client(profile=profile)
            logger.debug('__init__: len(client) = {}'.format(len(self.client)))
        except (IOError, IPython.parallel.TimeoutError):
            msg = 'Could not obtain an IPython parallel Client using profile "{}"'.format(profile)
            logger.error(msg)
            raise RuntimeError(msg)

    def __call__(self, task, args, message, batchsize):
        logger.debug('__call__: len(client) = {}'.format(len(self.client)))
        return list(parallel_map(self.client, task, args, message, batchsize))
