from multiprocessing import Pool
import logging
from functools import partial

from sage.all import ceil, floor


def binary_search(f, start, stop, param, predicate=lambda x, best: x <= best, *args, **kwds):
    """
    Searches for the best value in the interval [start,stop] depending on the given predicate.

    :param start: start of range to search
    :param stop: stop of range to search (exclusive)
    :param param: the parameter to modify when calling `f`
    :param predicate: comparison is performed by evaluating ``predicate(current, best)``
    """
    kwds[param] = stop
    D = {}
    D[stop] = f(*args, **kwds)
    best = D[stop]
    b = ceil((start + stop) / 2)
    direction = 0
    while True:
        if b not in D:
            kwds[param] = b
            D[b] = f(*args, **kwds)
        if b == start:
            best = D[start]
            break
        if not predicate(D[b], best):
            if direction == 0:
                start = b
                b = ceil((stop + b) / 2)
            else:
                stop = b
                b = floor((start + b) / 2)
        else:
            best = D[b]
            logging.getLogger("binsearch").debug("%s: %4d || %r" % (param, b, repr(best)))
            if b - 1 not in D:
                kwds[param] = b - 1
                D[b - 1] = f(*args, **kwds)
            if predicate(D[b - 1], best):
                stop = b
                b = floor((b + start) / 2)
                direction = 0
            else:
                if b + 1 not in D:
                    kwds[param] = b + 1
                    D[b + 1] = f(*args, **kwds)
                if not predicate(D[b + 1], best):
                    break
                else:
                    start = b
                    b = ceil((stop + b) / 2)
                    direction = 1

    for b in range(0, best[param])[::-1]:
        kwds[param] = b
        D[b] = f(*args, **kwds)
        if not predicate(D[b], best):
            break
        best = D[b]
        logging.getLogger("binsearch").debug("%s: %4d || %s" % (param, b, repr(best)))
    return best


def _batch_estimatef(f, x):
    y = f(x)
    logging.getLogger("batch").info(f"f: {f}")
    logging.getLogger("batch").info(f"x: {x}")
    logging.getLogger("batch").info(f"f(x): {repr(y)}")
    logging.getLogger("batch").info("")
    return y


def batch_estimate(params, algorithm, jobs=1, **kwds):
    from .lwe import LWEParameters

    if isinstance(params, LWEParameters):
        params = (params,)
    try:
        iter(algorithm)
    except TypeError:
        algorithm = (algorithm,)

    tasks = []

    for x in params:
        for f in algorithm:
            tasks.append((partial(f, **kwds), x))

    if jobs == 1:
        res = {}
        for f, x in tasks:
            y = _batch_estimatef(f, x)
            res[(f, x)] = y
    else:
        pool = Pool(jobs)
        res = pool.starmap(_batch_estimatef, tasks)
        res = dict([((f, x), res[i]) for i, (f, x) in enumerate(tasks)])

    return res
