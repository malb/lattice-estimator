from multiprocessing import Pool
from functools import partial

from sage.all import ceil, floor

from .io import Logging


class local_minimum:
    """
    An iterator context for finding a local minimum using binary search.

    We use the neighborhood of a point to decide the next direction to go into (gradient descent
    style), so the algorithm is not plain binary search (see ``update()`` function.)

    .. note :: We combine an iterator and a context to give the caller access to the result.
    """

    def __init__(
        self,
        start,
        stop,
        smallerf=lambda x, best: x <= best,
        log_level=5,
        suppress_bounds_warning=False,
    ):
        """
        Create a fresh local minimum search context.

        :param start: starting point
        :param stop:  end point (exclusive)
        :param smallerf: a function to decide if ``lhs`` is smaller than ``rhs``.
        :param suppress_bounds_warning: do not warn if a boundary is picked as optimal

        """
        self._suppress_bounds_warning = suppress_bounds_warning
        self._log_level = log_level
        self._start = start
        self._stop = stop - 1
        self._initial_bounds = (start, stop - 1)
        self._smallerf = smallerf
        # abs(self._direction) == 2: binary search step
        # abs(self._direction) == 1: gradient descent direction
        self._direction = -1  # going down
        self._last_x = None
        self._next_x = self._stop
        self._best = (None, None)

    def __enter__(self):
        """ """
        return self

    def __exit__(self, type, value, traceback):
        """ """
        pass

    def __iter__(self):
        """ """
        return self

    def __next__(self):
        if self._next_x is not None:
            self._last_x = self._next_x
            self._next_x = None
            return self._last_x
        else:
            if self._best[0] in self._initial_bounds and not self._suppress_bounds_warning:
                # We warn the user if the optimal solution is at the edge and thus possibly not optimal.
                Logging.log(
                    "bins",
                    self._log_level,
                    f'warning: "optimal" solution {self._best[0]} matches a bound ∈ {self._initial_bounds}.',
                )
            raise StopIteration

    @property
    def x(self):
        return self._best[0]

    @property
    def y(self):
        return self._best[1]

    def update(self, res):
        Logging.log("bins", self._log_level, f"({self._last_x}, {repr(res)})")

        # We got nothing yet
        if self._best[0] is None:
            self._best = self._last_x, res

        # We found something better
        if res is not False and self._smallerf(res, self._best[1]):
            # store it
            self._best = self._last_x, res

            # if it's a result of a long jump figure out the next direction
            if abs(self._direction) != 1:
                self._direction = -1
                self._next_x = self._last_x - 1
            # going down worked, so let's keep on doing that.
            elif self._direction == -1:
                self._direction = -2
                self._stop = self._last_x
                self._next_x = ceil((self._start + self._stop) / 2)
            # going up worked, so let's keep on doing that.
            elif self._direction == 1:
                self._direction = 2
                self._start = self._last_x
                self._next_x = floor((self._start + self._stop) / 2)
        else:
            # going downwards didn't help, let's try up
            if self._direction == -1:
                self._direction = 1
                self._next_x = self._last_x + 2
            # going up didn't help either, so we stop
            elif self._direction == 1:
                self._next_x = None
            # it got no better in a long jump, half the search space and try again
            elif self._direction == -2:
                self._start = self._last_x
                self._next_x = ceil((self._start + self._stop) / 2)
            elif self._direction == 2:
                self._stop = self._last_x
                self._next_x = floor((self._start + self._stop) / 2)

        # We are repeating ourselves, time to stop
        if self._next_x == self._last_x:
            self._next_x = None


def binary_search(
    f, start, stop, param, smallerf=lambda x, best: x <= best, log_level=5, *args, **kwds
):
    """
    Searches for the best value in the interval [start,stop] depending on the given comparison function.

    :param start: start of range to search
    :param stop: stop of range to search (exclusive)
    :param param: the parameter to modify when calling `f`
    :param smallerf: comparison is performed by evaluating ``smallerf(current, best)``
    """

    with local_minimum(start, stop + 1, smallerf=smallerf, log_level=log_level) as it:
        for x in it:
            kwds_ = dict(kwds)
            kwds_[param] = x
            it.update(f(*args, **kwds_))
        return it.y

    kwds[param] = stop
    D = {}
    D[stop] = f(*args, log_level=log_level + 1, **kwds)
    Logging.log("bins", log_level, f"{param}: {stop:4d} || {repr(D[stop])}")
    best = D[stop]
    b = ceil((start + stop) / 2)
    direction = 0
    while True:
        if b not in D:
            kwds[param] = b
            D[b] = f(*args, log_level=log_level + 1, **kwds)
            Logging.log("bins", log_level, f"{param}: {b:4d} || {repr(D[b])}")
        if b == start:
            best = D[start]
            break
        if not smallerf(D[b], best):
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
                D[b - 1] = f(*args, log_level=log_level + 1, **kwds)
                Logging.log("bins", log_level, f"{param}: {b-1:4d} || {repr(D[b-1])}")
            if smallerf(D[b - 1], best):
                best = D[b - 1]
                stop = b
                b = floor((b + start) / 2)
                direction = 0
            else:
                if b + 1 not in D:
                    kwds[param] = b + 1
                    D[b + 1] = f(*args, log_level=log_level + 1, **kwds)
                    Logging.log("bins", log_level, f"{param}: {b+1:4d} || {repr(D[b+1])}")
                if not smallerf(D[b + 1], best):
                    break
                else:
                    best = D[b + 1]
                    start = b
                    b = ceil((stop + b) / 2)
                    direction = 1
    return best


def _batch_estimatef(f, x, log_level=0, f_repr=None):
    y = f(x)
    if f_repr is None:
        f_repr = repr(f)
    Logging.log("batch", log_level, f"f: {f_repr}")
    Logging.log("batch", log_level, f"x: {x}")
    Logging.log("batch", log_level, f"f(x): {repr(y)}")
    return y


def f_name(f):
    try:
        return f.__name__
    except AttributeError:
        return repr(f)


def batch_estimate(params, algorithm, jobs=1, log_level=0, **kwds):
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
            tasks.append((partial(f, **kwds), x, log_level, f_name(f)))

    if jobs == 1:
        res = {}
        for f, x, lvl, f_repr in tasks:
            y = _batch_estimatef(f, x, lvl, f_repr)
            res[f_repr, x] = y
    else:
        pool = Pool(jobs)
        res = pool.starmap(_batch_estimatef, tasks)
        res = dict([((f_repr, x), res[i]) for i, (f, x, _, f_repr) in enumerate(tasks)])

    ret = dict()
    for f, x in res:
        ret[x] = ret.get(x, dict())
        ret[x][f] = res[f, x]

    return ret
