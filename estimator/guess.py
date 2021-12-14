# -*- coding: utf-8 -*-
"""
Generic multiplicative composition of guessing some components of the LWE secret and some LWE solving algorithm.

By "multiplicative" we mean that costs multiply rather than add. It is often possible to achieve
some form of additive composition, i.e. this strategy is rarely the most efficient.

"""
from sage.all import log, floor, ceil, binomial
from .util import local_minimum
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify


class guess_composition:
    def __init__(self, f):
        """
        Create a generic composition of guessing and `f`.
        """
        self.f = f
        self.__name__ = f"{f.__name__}+guessing"

    @classmethod
    def dense_solve(cls, f, params, log_level=5, **kwds):
        """
        Guess components of a dense secret then call `f`.

        :param f: Some object consuming `params` and outputting some `cost`
        :param params: LWE parameters.

        """
        base = params.Xs.bounds[1] - params.Xs.bounds[0] + 1

        baseline_cost = f(params, **kwds)

        max_zeta = floor(log(baseline_cost["rop"], base))

        with local_minimum(0, max_zeta, log_level=log_level) as it:
            for zeta in it:
                search_space = base ** zeta
                cost = f(params.updated(n=params.n - zeta), log_level=log_level + 1, **kwds)
                repeated_cost = cost.repeat(search_space)
                repeated_cost["zeta"] = zeta
                it.update(repeated_cost)
            return it.y

    @classmethod
    def gammaf(cls, n, h, zeta, base, g=lambda x: x):
        """
        Find optimal hamming weight for sparse guessing.

        Let `s` be a vector of dimension `n` where we expect `h` non-zero entries. We are ignoring `η-γ`
        components and are guessing `γ`. This succeeds with some probability given by ``prob_drop(n, h,
        ζ, γ)``. Exhaustively searching the guesses takes `binomial(n, γ) ⋅ b^γ` steps where `b` is the
        number of non-zero values in a component of `s`. We call a `γ` optimal if it minimizes the
        overall number of repetitions that need to be performed to succeed with probability 99%.

        :param n: vector dimension
        :param h: hamming weight of the vector
        :param zeta: number of ignored + guesses components
        :param base: number of possible non-zero scalars
        :param g: We do not consider search space directly by `g()` applied to it (think time-memory
                  trade-offs).
        :returns: (number of repetitions, γ, size of the search space, probability of success)

        """
        if not zeta:
            return 1, 0, 0, 1.0

        search_space = 0
        gamma = 0
        probability = 0
        best = None, None, None, None
        while gamma < min(h, zeta):
            probability += prob_drop(n, h, zeta, fail=gamma)
            search_space += binomial(zeta, gamma) * base ** gamma
            repeat = prob_amplify(0.99, probability) * g(search_space)
            if best[0] is None or repeat < best[0]:
                best = repeat, gamma, search_space, probability
                gamma += 1
            else:
                break
        return best

    @classmethod
    def sparse_solve(cls, f, params, log_level=5, **kwds):
        """
        Guess components of a sparse secret then call `f`.

        :param f: Some object consuming `params` and outputting some `cost`
        :param params: LWE parameters.
        """
        base = params.Xs.bounds[1] - params.Xs.bounds[0]  # we exclude zero
        h = ceil(len(params.Xs) * params.Xs.density)  # nr of non-zero entries

        with local_minimum(0, params.n - 40, log_level=log_level) as it:
            for zeta in it:
                single_cost = f(params.updated(n=params.n - zeta), log_level=log_level + 1, **kwds)
                repeat, gamma, search_space, probability = cls.gammaf(params.n, h, zeta, base)
                cost = single_cost.repeat(repeat)
                cost["zeta"] = zeta
                cost["|S|"] = search_space
                cost["prop"] = probability
                it.update(cost)
            return it.y

    def __call__(self, params, log_level=5, **kwds):
        """
        Guess components of a secret then call `f`.

        :param params: LWE parameters.

        EXAMPLE::

            >>> from estimator import *
            >>> from estimator.guess import guess_composition
            >>> guess_composition(primal_usvp)(Kyber512.updated(Xs=ND.SparseTernary(512, 16)))
            rop: ≈2^102.0, red: ≈2^102.0, δ: 1.008079, β: 130, d: 474, tag: usvp, ↻: ≈2^31.8, ζ: 235, |S|: 1, ...

        Compare::

            >>> primal_hybrid(Kyber512.updated(Xs=ND.SparseTernary(512, 16)))
            rop: ≈2^85.8, red: ≈2^84.9, svp: ≈2^84.7, β: 104, η: 2, ζ: 368, |S|: ≈2^91.0, d: 311, prob: ≈2^-20...

        """
        if params.Xs.is_sparse:
            return self.sparse_solve(self.f, params, log_level, **kwds)
        else:
            return self.dense_solve(self.f, params, log_level, **kwds)
