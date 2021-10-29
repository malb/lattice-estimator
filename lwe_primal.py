# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE.
"""
from functools import partial
from multiprocessing import Pool

from sage.all import oo, ceil, sqrt, log, RR
from .reduction import BKZ
from .util import binary_search
from .lwe import LWEParameters


class PrimalUSVP:
    @staticmethod
    def _scale_factor(Xs, Xe):
        scale = RR(1)
        if Xs < Xe:
            scale = Xe.stddev / Xs.stddev
        return scale

    @staticmethod
    def _solve_for_d(params, m, beta, kannan_coeff, scale):
        """
        Find smallest d ∈ [n,m] to satisfy uSVP condition.

        If no such d exists, return the upper bound m.
        """
        # Find the smallest d ∈ [n,m] s.t. a*d^2 + b*d + c >= 0
        delta = BKZ.delta(beta)
        a = -log(delta)
        C = log(params.Xe.stddev ** 2 * (beta - 1) + kannan_coeff ** 2) / 2.0
        b = log(delta) * (2 * beta - 1) + log(params.q) - C
        c = log(kannan_coeff) + params.n * log(scale) - (params.n + 1) * log(params.q)
        n = params.n
        if a * n * n + b * n + c >= 0:  # trivial case
            return n

        # solve for ad^2 + bd + c == 0
        disc = b * b - 4 * a * c  # the discriminant
        if disc < 0:  # no solution, return m
            return m

        # compute the two solutions
        d1 = (-b + sqrt(disc)) / (2 * a)
        d2 = (-b - sqrt(disc)) / (2 * a)
        if a > 0:  # the only possible solution is ceiling(d2)
            return min(m, ceil(d2))

        # the case a<=0:
        # if n is to the left of d1 then the first solution is ceil(d1)
        if n <= d1:
            return min(m, ceil(d1))

        # otherwise, n must be larger than d2 (since an^2+bn+c<0) so no solution
        return m

    @staticmethod
    def cost_gsa(
        beta: int,
        params: LWEParameters,
        m: int = oo,
        kannan_coeff=None,
        d=None,
        reduction_cost_model=BKZ.default,
    ):

        delta = BKZ.delta(beta)
        scale = PrimalUSVP._scale_factor(params.Xs, params.Xe)
        m = min(2 * ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        kannan_coeff = params.Xe.stddev if kannan_coeff is None else kannan_coeff
        d = PrimalUSVP._solve_for_d(params, m, beta, kannan_coeff, scale) if d is None else d
        assert d <= m + 1

        lhs = log(sqrt(params.Xe.stddev ** 2 * (beta - 1) + kannan_coeff ** 2))
        rhs = RR(
            log(delta) * (2 * beta - d - 1)
            + (log(kannan_coeff) + log(scale) * params.n + log(params.q) * (d - params.n - 1)) / d
        )

        return BKZ.cost(reduction_cost_model, beta, d, predicate=lhs <= rhs)

    @staticmethod
    def cost_simulator(
        beta: int,
        params: LWEParameters,
        simulator,
        m: int = oo,
        kannan_coeff=None,
        d=None,
        reduction_cost_model=BKZ.default,
    ):
        delta = BKZ.delta(beta)
        m = min(2 * ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        d = m + 1 if d is None else d
        scale = PrimalUSVP._scale_factor(params.Xe, params.Xs)
        kannan_coeff = params.Xe.stddev if kannan_coeff is None else kannan_coeff

        def f(d):
            r = simulator(d, params.n, params.q, beta=beta, scale=scale, kannan_coeff=kannan_coeff)
            lhs = params.Xe.stddev ** 2 * (beta - 1) + kannan_coeff ** 2
            if r[d - beta] > lhs:
                cost = BKZ.cost(reduction_cost_model, beta, d)
            else:
                cost = BKZ.cost(reduction_cost_model, beta, d, predicate=False)
            return cost

        if f(d)["red"] < oo:
            cost = binary_search(
                f,
                start=params.n,
                stop=d,
                param="d",
                predicate=lambda x, best: x["red"] <= best["red"],
            )
            return cost
        else:
            return f(d)

    def __call__(
        self,
        params: LWEParameters,
        success_probability=0.99,
        kannan_coeff=None,
        d=None,
        reduction_cost_model=BKZ.default,
        bkz_model="gsa",
        **kwds,
    ):
        params = LWEParameters.normalize(params)

        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

        if bkz_model == "gsa":
            f = self.cost_gsa
            start = 40
            stop = 2 * params.n
        elif str(bkz_model).lower() != "gsa":
            cost_gsa = self(
                params,
                kannan_coeff=kannan_coeff,
                reduction_cost_model=reduction_cost_model,
                bkz_model="gsa",
                **kwds,
            )
            start = cost_gsa["beta"] - 32
            stop = cost_gsa["beta"] + 128
            f = partial(self.cost_simulator, simulator=bkz_model)
        else:
            raise NotImplementedError(f"{bkz_model} unknown")

        cost = binary_search(
            f,
            start=start,
            stop=stop,
            param="beta",
            predicate=lambda x, best: x["red"] <= best["red"],
            params=params,
            kannan_coeff=kannan_coeff,
            d=d,
            m=m,
            reduction_cost_model=reduction_cost_model,
        )

        return cost

    def __repr__(self):
        return "primal_usvp"


primal_usvp = PrimalUSVP()
primal_usvp_cn11 = partial(primal_usvp, bkz_model="cn11")


def _fpp(f, x):
    y = f(x)
    print(f)
    print(x)
    print(y)
    print("")
    return y


def batch_estimate(params, algorithm, jobs=1, **kwds):
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
            y = _fpp(x)
            res[(f, x)] = y
    else:
        pool = Pool(jobs)
        res = pool.starmap(_fpp, tasks)
        res = dict([((f, x), res[i]) for i, (f, x) in enumerate(tasks)])

    return res
