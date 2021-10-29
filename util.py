from sage.all import ceil, floor
import logging


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
            logging.getLogger("binsearch").debug("%s: %4d || %s" % (param, b, best))
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

    for b in range(40, best[param])[::-1]:
        kwds[param] = b
        D[b] = f(*args, **kwds)
        if not predicate(D[b], best):
            break
        best = D[b]
    return best
