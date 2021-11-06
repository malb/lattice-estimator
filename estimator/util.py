from multiprocessing import Pool
from functools import partial

from sage.all import ceil, floor

from .io import Logging


def binary_search(
    f, start, stop, param, predicate=lambda x, best: x <= best, log_level=5, *args, **kwds
):
    """
    Searches for the best value in the interval [start,stop] depending on the given predicate.

    :param start: start of range to search
    :param stop: stop of range to search (exclusive)
    :param param: the parameter to modify when calling `f`
    :param predicate: comparison is performed by evaluating ``predicate(current, best)``
    """
    bounds = (start, stop)

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
            Logging.log("bins", log_level, f"{param}: {b:4d} || {repr(best)}")
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

    for b in range(bounds[0], best[param])[::-1]:
        kwds[param] = b
        D[b] = f(*args, **kwds)
        if not predicate(D[b], best):
            break
        best = D[b]
        Logging.log("bins", log_level, f"{param}: {b:4d} || {repr(best)}")
    return best


def _batch_estimatef(f, x, log_level=0, f_repr=None):
    y = f(x)
    if f_repr is None:
        f_repr = repr(f)
    Logging.log("batch", log_level, f"f: {f_repr}")
    Logging.log("batch", log_level, f"x: {x}")
    Logging.log("batch", log_level, f"f(x): {repr(y)}")
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
            try:
                f_repr = f.__name__
            except AttributeError:
                f_repr = repr(f)
            tasks.append((partial(f, **kwds), x, 0, f_repr))

    if jobs == 1:
        res = {}
        for f, x, lvl, f_repr in tasks:
            y = _batch_estimatef(f, x, lvl, f_repr)
            res[(f_repr, x)] = y
    else:
        pool = Pool(jobs)
        res = pool.starmap(_batch_estimatef, tasks)
        res = dict([((f_repr, x), res[i]) for i, (f, x, _, f_repr) in enumerate(tasks)])

    ret = dict()
    for f, x in res:
        ret[x] = ret.get(x, dict())
        ret[x][f] = res[f, x]

    return ret
