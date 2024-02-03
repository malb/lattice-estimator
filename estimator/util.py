import itertools as it
from multiprocessing import Pool
from functools import partial
from dataclasses import dataclass, field
from typing import Any, Callable, NamedTuple

from sage.all import ceil, floor, log, oo, RR, cached_function, zeta

from .io import Logging
from .lwe_parameters import LWEParameters
from .sis_parameters import SISParameters
from .conf import max_n_cache


def log2(x):
    return log(x, 2.0)


@cached_function
def zeta_prime(x):
    h = 1e-5
    return RR((zeta(x+h) - zeta(x-h)))/(2*h)


# Low beta Gaussian Heuristic constant for use in NTRU Dense sublattice estimation.
gh_constant = {1: 0.00000, 2: -0.50511, 3: -0.46488, 4: -0.39100, 5: -0.29759, 6: -0.24880, 7: -0.21970, 8: -0.15748,
               9: -0.14673, 10: -0.07541, 11: -0.04870, 12: -0.01045, 13: 0.02298, 14: 0.04212, 15: 0.07014,
               16: 0.09205, 17: 0.12004, 18: 0.14988, 19: 0.17351, 20: 0.18659, 21: 0.20971, 22: 0.22728, 23: 0.24951,
               24: 0.26313, 25: 0.27662, 26: 0.29430, 27: 0.31399, 28: 0.32494, 29: 0.34796, 30: 0.36118, 31: 0.37531,
               32: 0.39056, 33: 0.39958, 34: 0.41473, 35: 0.42560, 36: 0.44222, 37: 0.45396, 38: 0.46275, 39: 0.47550,
               40: 0.48889, 41: 0.50009, 42: 0.51312, 43: 0.52463, 44: 0.52903, 45: 0.53930, 46: 0.55289, 47: 0.56343,
               48: 0.57204, 49: 0.58184, 50: 0.58852}

# Low beta \alpha_\beta quantity as defined in [AC:DucWoe21] for use in NTRU Dense subblattice estimation.
small_slope_t8 = {2: 0.04473, 3: 0.04472, 4: 0.04402, 5: 0.04407, 6: 0.04334, 7: 0.04326, 8: 0.04218, 9: 0.04237,
                  10: 0.04144, 11: 0.04054, 12: 0.03961, 13: 0.03862, 14: 0.03745, 15: 0.03673, 16: 0.03585,
                  17: 0.03477, 18: 0.03378, 19: 0.03298, 20: 0.03222, 21: 0.03155, 22: 0.03088, 23: 0.03029,
                  24: 0.02999, 25: 0.02954, 26: 0.02922, 27: 0.02891, 28: 0.02878, 29: 0.02850, 30: 0.02827,
                  31: 0.02801, 32: 0.02786, 33: 0.02761, 34: 0.02768, 35: 0.02744, 36: 0.02728, 37: 0.02713,
                  38: 0.02689, 39: 0.02678, 40: 0.02671, 41: 0.02647, 42: 0.02634, 43: 0.02614, 44: 0.02595,
                  45: 0.02583, 46: 0.02559, 47: 0.02534, 48: 0.02514, 49: 0.02506, 50: 0.02493, 51: 0.02475,
                  52: 0.02454, 53: 0.02441, 54: 0.02427, 55: 0.02407, 56: 0.02393, 57: 0.02371, 58: 0.02366,
                  59: 0.02341, 60: 0.02332}


@dataclass
class LazyEvaluation:
    f: Callable
    max_n_cache: int
    eval: list = field(default_factory=lambda: [])

    def __getitem__(self, key):
        if not self.eval:
            self.eval = [self.f(i) for i in range(self.max_n_cache + 1)]

        return self.eval[key]


zeta_precomputed = LazyEvaluation(lambda i: RR(zeta(i)) if i != 1 else RR(oo), max_n_cache)
zeta_prime_precomputed = LazyEvaluation(zeta_prime, max_n_cache)


class Bounds(NamedTuple):
    low: Any
    high: Any


class local_minimum_base:
    """
    An iterator context for finding a local minimum using binary search.

    We use the immediate neighborhood of a point to decide the next direction to go into (gradient
    descent style), so the algorithm is not plain binary search (see ``update()`` function.)

    .. note :: We combine an iterator and a context to give the caller access to the result.
    """

    def __init__(
        self,
        start,
        stop,
        smallerf=lambda x, best: x <= best,
        suppress_bounds_warning=False,
        log_level=5,
    ):
        """
        Create a fresh local minimum search context.

        :param start: starting point
        :param stop:  end point (exclusive)
        :param smallerf: a function to decide if ``lhs`` is smaller than ``rhs``.
        :param suppress_bounds_warning: do not warn if a boundary is picked as optimal

        """

        if stop < start:
            raise ValueError(f"Incorrect bounds {start} > {stop}.")

        self._suppress_bounds_warning = suppress_bounds_warning
        self._log_level = log_level
        self._start = start
        self._stop = stop - 1
        self._initial_bounds = Bounds(start, stop - 1)
        self._smallerf = smallerf
        # abs(self._direction) == 2: binary search step
        # abs(self._direction) == 1: gradient descent direction
        self._direction = -1  # going down
        self._last_x = None
        self._next_x = self._stop
        self._best = Bounds(None, None)
        self._all_x = set()

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

        if (
            self._next_x is not None
            and self._next_x not in self._all_x
            and self._initial_bounds.low <= self._next_x <= self._initial_bounds.high
        ):
            # we've not been told to abort
            # we're not looping
            # we're in bounds
            self._last_x = self._next_x
            self._next_x = None
            return self._last_x

        if self._best.low in self._initial_bounds and not self._suppress_bounds_warning:
            # We warn the user if the optimal solution is at the edge and thus possibly not optimal.
            msg = (
                f'warning: "optimal" solution {self._best.low} matches a bound ∈ {self._initial_bounds}.',
            )
            Logging.log("bins", self._log_level, msg)

        raise StopIteration

    @property
    def x(self):
        return self._best.low

    @property
    def y(self):
        return self._best.high

    def update(self, res):
        """

        TESTS:

        We keep cache old inputs in ``_all_x`` to prevent infinite loops::

            >>> from estimator.util import binary_search
            >>> from estimator.cost import Cost
            >>> f = lambda x, log_level=1: Cost(rop=1) if x >= 19 else Cost(rop=2)
            >>> binary_search(f, 10, 30, "x")
            rop: 1

        """

        Logging.log("bins", self._log_level, f"({self._last_x}, {repr(res)})")

        self._all_x.add(self._last_x)

        # We got nothing yet
        if self._best.low is None:
            self._best = Bounds(self._last_x, res)

        # We found something better
        if res is not False and self._smallerf(res, self._best.high):
            # store it
            self._best = Bounds(self._last_x, res)

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


class local_minimum(local_minimum_base):
    """
    An iterator context for finding a local minimum using binary search.

    We use the neighborhood of a point to decide the next direction to go into (gradient descent
    style), so the algorithm is not plain binary search (see ``update()`` function.)

    We also zoom out by a factor ``precision``, find an approximate local minimum and then
    search the neighbourhood for the smallest value.

    .. note :: We combine an iterator and a context to give the caller access to the result.

    """

    def __init__(
        self,
        start,
        stop,
        precision=1,
        smallerf=lambda x, best: x <= best,
        suppress_bounds_warning=False,
        log_level=5,
    ):
        """
        Create a fresh local minimum search context.

        :param start: starting point
        :param stop:  end point (exclusive)
        :param precision: only consider every ``precision``-th value in the main loop
        :param smallerf: a function to decide if ``lhs`` is smaller than ``rhs``.
        :param suppress_bounds_warning: do not warn if a boundary is picked as optimal

        """
        self._precision = precision
        self._orig_bounds = (start, stop)
        start = ceil(start / precision)
        stop = floor(stop / precision)
        local_minimum_base.__init__(self, start, stop, smallerf, suppress_bounds_warning, log_level)

    def __next__(self):
        x = local_minimum_base.__next__(self)
        return x * self._precision

    @property
    def x(self):
        return self._best.low * self._precision

    @property
    def neighborhood(self):
        """
        An iterator over the neighborhood of the currently best value.
        """

        start_bound, stop_bound = self._orig_bounds
        start = max(start_bound, self.x - self._precision)
        stop = min(stop_bound, self.x + self._precision)
        return range(start, stop)


class early_abort_range:
    """
    An iterator context for finding a local minimum using linear search.

    .. note :: We combine an iterator and a context to give the caller access to the result.
    """

    # TODO: unify whether we like contexts or not

    def __init__(
        self,
        start,
        stop=oo,
        step=1,
        smallerf=lambda x, best: x <= best,
        suppress_bounds_warning=False,
        log_level=5,
    ):
        """
        Create a fresh local minimum search context.

        :param start: starting point
        :param stop:  end point (exclusive, optional)
        :param step:  step size
        :param smallerf: a function to decide if ``lhs`` is smaller than ``rhs``.
        :param suppress_bounds_warning: do not warn if a boundary is picked as optimal

        """

        if stop < start:
            raise ValueError(f"Incorrect bounds {start} > {stop}.")

        self._suppress_bounds_warning = suppress_bounds_warning
        self._log_level = log_level
        self._start = start
        self._step = step
        self._stop = stop
        self._smallerf = smallerf
        self._last_x = None
        self._next_x = self._start
        self._best = Bounds(None, None)

    def __iter__(self):
        """ """
        return self

    def __next__(self):
        if self._next_x is None:
            raise StopIteration
        if self._next_x >= self._stop:
            raise StopIteration

        self._last_x = self._next_x
        self._next_x += self._step
        return self._last_x, self

    @property
    def x(self):
        return self._best.low

    @property
    def y(self):
        return self._best.high

    def update(self, res):
        """ """
        Logging.log("lins", self._log_level, f"({self._last_x}, {repr(res)})")

        if self._best.low is None:
            self._best = Bounds(self._last_x, res)
            return

        if res is False:
            self._next_x = None
        elif self._smallerf(res, self._best.high):
            self._best = Bounds(self._last_x, res)
        else:
            self._next_x = None


def binary_search(
    f, start, stop, param, step=1, smallerf=lambda x, best: x <= best, log_level=5, *args, **kwds
):
    """
    Searches for the best value in the interval [start,stop] depending on the given comparison function.

    :param start: start of range to search (inclusive)
    :param stop: stop of range to search (inclusive)
    :param param: the parameter to modify when calling `f`
    :param smallerf: comparison is performed by evaluating ``smallerf(current, best)``
    :param step: initially only consider every `step`-th value
    """

    with local_minimum(start, stop + 1, step, smallerf=smallerf, log_level=log_level) as it:
        for x in it:
            kwds_ = dict(kwds)
            kwds_[param] = x
            it.update(f(*args, **kwds_))

        for x in it.neighborhood:
            kwds_ = dict(kwds)
            kwds_[param] = x
            it.update(f(*args, **kwds_))

        return it.y


def _batch_estimatef(f, x, log_level=0, f_repr=None, catch_exceptions=True):
    try:
        y = f(x)
    except Exception as e:
        if catch_exceptions:
            print(f"Algorithm {f_repr} on {x} failed with {e}")
            return None
        else:
            raise e

    if f_repr is None:
        f_repr = repr(f)

    Logging.log("batch", log_level, f"f: {f_repr}")
    Logging.log("batch", log_level, f"x: {x}")
    Logging.log("batch", log_level, f"f(x): {y!r}")

    return y


def f_name(f):
    try:
        return f.__name__
    except AttributeError:
        return repr(f)


class Task(NamedTuple):
    f: Callable
    x: LWEParameters
    log_level: int
    f_name: str
    catch_exceptions: bool


@dataclass(frozen=True)
class TaskResults:
    _map: dict

    def __getitem__(self, params):
        return {
            task.f_name: result
            for task, result in self._map.items()
            if task.x == params and result is not None
        }


def batch_estimate(params, algorithm, jobs=1, log_level=0, catch_exceptions=True, **kwds):
    """
    Run estimates for all algorithms for all parameters.

    :param params: (List of) LWE parameters.
    :param algorithm: (List of) algorithms.
    :param jobs: Use multiple threads in parallel.
    :param log_level:
    :param catch_exceptions: When an estimate fails, just print a warning.

    Example::

        >>> from estimator import LWE
        >>> from estimator.schemes import Kyber512
        >>> _ = batch_estimate(Kyber512, [LWE.primal_usvp, LWE.primal_bdd])
        >>> _ = batch_estimate(Kyber512, [LWE.primal_usvp, LWE.primal_bdd], jobs=2)

    """

    if isinstance(params, LWEParameters) or isinstance(params, SISParameters):
        params = (params,)
    if not hasattr(algorithm, "__iter__"):
        algorithm = (algorithm,)
    tasks = [
        Task(partial(f, **kwds), x, log_level, f_name(f), catch_exceptions)
        for f, x in it.product(algorithm, params)
    ]

    if jobs == 1:
        results = [_batch_estimatef(*task) for task in tasks]
    else:
        with Pool(jobs) as pool:
            results = pool.starmap(_batch_estimatef, tasks)

    return TaskResults(dict(zip(tasks, results)))
