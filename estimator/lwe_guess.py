# -*- coding: utf-8 -*-
"""
Generic multiplicative composition of guessing some components of the LWE secret and some LWE solving algorithm.

By "multiplicative" we mean that costs multiply rather than add. It is often possible to achieve
some form of additive composition, i.e. this strategy is rarely the most efficient.

"""

from sage.all import binomial, ceil, e, exp, floor, log, oo, pi, round, RR, sqrt, ZZ

from .conf import mitm_opt
from .cost import Cost
from .errors import InsufficientSamplesError, OutOfBoundsError
from .lwe_parameters import LWEParameters
from .prob import amplify as prob_amplify, drop as prob_drop, amplify_sigma
from .util import local_minimum, log2
from .nd import sigmaf, SparseTernary


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

        if baseline_cost["rop"] == oo:
            # yeah, no
            return baseline_cost

        max_zeta = min(floor(log(baseline_cost["rop"], base)), params.n)

        with local_minimum(0, max_zeta, log_level=log_level) as it:
            for zeta in it:
                search_space = base**zeta
                cost = f(params.updated(n=params.n - zeta), log_level=log_level + 1, **kwds)
                if cost["rop"] == oo:
                    return Cost(rop=oo)
                repeated_cost = cost.repeat(search_space)
                repeated_cost["zeta"] = zeta
                it.update(repeated_cost)
            return it.y if it.y else Cost(rop=oo)

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
        if h == 0 or not zeta:
            return 1, 0, 0, 1.0

        search_space = 0
        gamma = 0
        probability = 0
        best = None, None, None, None
        while gamma < min(h, zeta):
            probability += prob_drop(n, h, zeta, fail=gamma)
            search_space += binomial(zeta, gamma) * base**gamma
            repeat = prob_amplify(0.99, probability) * g(search_space)
            if best[0] is not None and repeat >= best[0]:
                break
            best = repeat, gamma, search_space, probability
            gamma += 1
        return best

    @classmethod
    def sparse_solve(cls, f, params, log_level=5, **kwds):
        """
        Guess components of a sparse secret then call `f`.

        :param f: Some object consuming `params` and outputting some `cost`
        :param params: LWE parameters.
        """
        base = params.Xs.bounds[1] - params.Xs.bounds[0]  # we exclude zero
        h = params.Xs.hamming_weight

        with local_minimum(0, params.n - 40, log_level=log_level) as it:
            for zeta in it:
                single_cost = f(params.updated(n=params.n - zeta), log_level=log_level + 1, **kwds)
                if single_cost["rop"] == oo:
                    return Cost(rop=oo)
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
            >>> from estimator.lwe_guess import guess_composition
            >>> guess_composition(LWE.primal_usvp)(schemes.Kyber512.updated(Xs=ND.SparseTernary(16)))
            rop: ≈2^102.1, red: ≈2^102.1, δ: 1.008011, β: 132, d: 461, tag: usvp, ↻: ≈2^34.9, ζ: 252, |S|: 1, ...

        Compare::

            >>> LWE.primal_hybrid(schemes.Kyber512.updated(Xs=ND.SparseTernary(16)))
            rop: ≈2^85.8, red: ≈2^84.8, svp: ≈2^84.8, β: 105, η: 2, ζ: 364, |S|: ≈2^85.0, d: 317, prob: ≈2^-23.4, ↻:...

        """
        params = LWEParameters.normalize(params)
        solve = self.sparse_solve if params.Xs.is_sparse else self.dense_solve
        return solve(self.f, params, log_level, **kwds)


class ExhaustiveSearch:
    def __call__(self, params: LWEParameters, success_probability=0.99, quantum: bool = False):
        """
        Estimate cost of solving LWE via exhaustive search.

        :param params: LWE parameters
        :param success_probability: the targeted success probability
        :param quantum: use estimate for quantum computer (we simply take the square root of the search space)
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: memory requirement in integers mod q.
        - ``m``: Required number of samples to distinguish the correct solution with high probability.

        EXAMPLE::

            >>> from estimator import *
            >>> from estimator.lwe_guess import exhaustive_search
            >>> params = LWE.Parameters(n=64, q=2**40, Xs=ND.Binary, Xe=ND.DiscreteGaussian(3.2))
            >>> exhaustive_search(params)
            rop: ≈2^73.6, mem: ≈2^72.6, m: 397.198
            >>> params = LWE.Parameters(n=1024, q=2**40, Xs=ND.SparseTernary(32), Xe=ND.DiscreteGaussian(3.2))
            >>> exhaustive_search(params)
            rop: ≈2^413.9, mem: ≈2^412.9, m: ≈2^11.1

        """
        params = LWEParameters.normalize(params)

        # there are two stages: enumeration and distinguishing, so we split up the success_probability
        probability = sqrt(success_probability)

        try:
            size = params.Xs.support_size(probability)
        except NotImplementedError:
            # not achieving required probability with search space
            # given our settings that means the search space is huge
            # so we approximate the cost with oo
            return Cost(rop=oo, mem=oo, m=1)

        if quantum is True:
            size = size.sqrt()

        # set m according to [ia.cr/2020/515]
        sigma = params.Xe.stddev / params.q
        m_required = RR(
            8 * exp(4 * pi * pi * sigma * sigma) * (log(size) - log(log(1 / probability)))
        )

        if params.m < m_required:
            raise InsufficientSamplesError(
                f"Exhaustive search: Need {m_required} samples but only {params.m} available."
            )

        # we can compute A*s for all candidate s in time 2*size*m using
        # (the generalization [ia.cr/2021/152] of) the recursive algorithm
        # from [ia.cr/2020/515]
        cost = 2 * size * m_required

        ret = Cost(rop=cost, mem=cost / 2, m=m_required)
        return ret.sanity_check()

    __name__ = "exhaustive_search"


exhaustive_search = ExhaustiveSearch()


class MITM:

    locality = 0.05

    def X_range(self, nd):
        if nd.is_bounded:
            a, b = nd.bounds
            return b - a + 1, 1.0
        else:
            # setting fraction=0 to ensure that support size does not
            # throw error. we'll take the probability into account later
            rng = nd.resize(1).support_size(0.0)
            return rng, nd.gaussian_tail_prob

    def local_range(self, center):
        return ZZ(floor((1 - self.locality) * center)), ZZ(ceil((1 + self.locality) * center))

    def mitm_analytical(self, params: LWEParameters, success_probability=0.99):
        nd_rng, nd_p = self.X_range(params.Xe)
        delta = nd_rng / params.q  # possible error range scaled

        sd_rng, sd_p = self.X_range(params.Xs)

        # determine the number of elements in the tables depending on splitting dim
        n = params.n
        k = round(n / (2 - delta))
        # we could now call self.cost with this k, but using our model below seems
        # about 3x faster and reasonably accurate

        if params.Xs.is_sparse:
            h = params.Xs.hamming_weight
            if type(params.Xs) is SparseTernary:
                # split optimally and compute the probability of this event
                success_probability_ = params.Xs.split_probability(k)
            else:
                split_h = (h * k / n).round("down")
                # Assume each coefficient is sampled i.i.d.:
                success_probability_ = (
                    binomial(k, split_h) * binomial(n - k, h - split_h) / binomial(n, h)
                )

            logT = RR(h * (log2(n) - log2(h) + log2(sd_rng - 1) + log2(e))) / (2 - delta)
            logT -= RR(log2(h) / 2)
            logT -= RR(h * h * log2(e) / (2 * n * (2 - delta) ** 2))
        else:
            success_probability_ = 1.0
            logT = k * log(sd_rng, 2)

        m_required = max(1, round(logT + log2(logT)))
        if params.m < m_required:
            raise InsufficientSamplesError(
                f"MITM: Need {m_required} samples but only {params.m} available."
            )

        # since m = logT + loglogT and rop = T*m, we have rop=2^m
        ret = Cost(rop=RR(2**m_required), mem=2**logT * m_required, m=m_required, k=ZZ(k))
        repeat = prob_amplify(
            success_probability, sd_p**n * nd_p**m_required * success_probability_
        )
        return ret.repeat(times=repeat)

    def cost(
        self,
        params: LWEParameters,
        k: int,
        success_probability=0.99,
    ):
        nd_rng, nd_p = self.X_range(params.Xe)
        delta = nd_rng / params.q  # possible error range scaled

        sd_rng, sd_p = self.X_range(params.Xs)
        n = params.n

        if params.Xs.is_sparse:
            h = params.Xs.hamming_weight
            # we assume the hamming weight to be distributed evenly across the two parts
            # if not we can rerandomize on the coordinates and try again -> repeat
            if type(params.Xs) is SparseTernary:
                sec_tab, sec_sea = params.Xs.split_balanced(k)
                size_tab = sec_tab.support_size()
                size_sea = sec_sea.support_size()
            else:
                # Assume each coefficient is sampled i.i.d.:
                split_h = (h * k / n).round("down")
                size_tab = RR((sd_rng - 1) ** split_h * binomial(k, split_h))
                size_sea = RR((sd_rng - 1) ** (h - split_h) * binomial(n - k, h - split_h))

            success_probability_ = size_tab * size_sea / params.Xs.support_size()
        else:
            size_tab = sd_rng**k
            size_sea = sd_rng ** (n - k)
            success_probability_ = 1

        # we set m such that it approximately minimizes the search cost per query as
        # a reasonable starting point and then optimize around it
        m_ = ceil(max(log2(size_tab) + log2(log2(size_tab)), 1))
        a, b = self.local_range(m_)
        with local_minimum(a, b, smallerf=lambda x, best: x[1] <= best[1]) as it:
            for m in it:
                # for search we effectively build a second table and for each entry, we expect
                # 2^( m * 4 * B / q) = 2^(delta * m) table look ups + a l_oo computation (costing m)
                # for every hit in the table (which has probability T/2^m)
                cost = (m, size_sea * (2 * m + 2 ** (delta * m) * (1 + size_tab * m / 2**m)))
                it.update(cost)
            m, cost_search = it.y
        m = min(m, params.m)

        # building the table costs 2*T*m using the generalization [ia.cr/2021/152] of
        # the recursive algorithm from [ia.cr/2020/515]
        cost_table = size_tab * 2 * m

        ret = Cost(rop=(cost_table + cost_search), m=m, k=k)
        ret["mem"] = size_tab * (k + m) + size_sea * (n - k + m)
        repeat = prob_amplify(success_probability, sd_p**n * nd_p**m * success_probability_)
        return ret.repeat(times=repeat)

    def __call__(self, params: LWEParameters, success_probability=0.99, optimization=mitm_opt):
        """
        Estimate cost of solving LWE via Meet-In-The-Middle attack.

        :param params: LWE parameters
        :param success_probability: the targeted success probability
        :param model: Either "analytical" (faster, default) or "numerical" (more accurate)
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: memory requirement in integers mod q.
        - ``m``: Required number of samples to distinguish the correct solution with high probability.
        - ``k``: Splitting dimension.
        - ``↻``: Repetitions required to achieve targeted success probability

        EXAMPLE::

            >>> from estimator import *
            >>> from estimator.lwe_guess import mitm
            >>> params = LWE.Parameters(n=64, q=2**40, Xs=ND.Binary, Xe=ND.DiscreteGaussian(3.2))
            >>> mitm(params)
            rop: ≈2^37.0, mem: ≈2^37.2, m: 37, k: 32, ↻: 1
            >>> mitm(params, optimization="numerical")
            rop: ≈2^39.2, m: 36, k: 32, mem: ≈2^39.1, ↻: 1
            >>> params = LWE.Parameters(n=1024, q=2**40, Xs=ND.SparseTernary(32), Xe=ND.DiscreteGaussian(3.2))
            >>> mitm(params)
            rop: ≈2^217.8, mem: ≈2^210.2, m: ≈2^15.5, k: 512, ↻: 226
            >>> mitm(params, optimization="numerical")
            rop: ≈2^215.6, m: ≈2^15.5, k: 512, mem: ≈2^208.6, ↻: 226

        """
        Cost.register_impermanent(rop=True, mem=False, m=True, k=False)

        params = LWEParameters.normalize(params)

        nd_rng = self.X_range(params.Xe)[0]
        if nd_rng >= params.q:
            # MITM attacks cannot handle an error this large.
            return Cost(rop=oo, mem=oo, m=0, k=0)

        if "analytical" in optimization:
            return self.mitm_analytical(params=params, success_probability=success_probability)
        elif "numerical" in optimization:
            with local_minimum(1, params.n - 1) as it:
                for k in it:
                    cost = self.cost(k=k, params=params, success_probability=success_probability)
                    it.update(cost)
                ret = it.y
                # if the noise is large, the curve might not be convex, so the above minimum
                # is not correct. Interestingly, in these cases, it seems that k=1 might be smallest
                ret1 = self.cost(k=1, params=params, success_probability=success_probability)
                return min(ret, ret1)
        else:
            raise ValueError("Unknown optimization method for MITM.")

    __name__ = "mitm"


mitm = MITM()


class Distinguisher:
    def __call__(self, params: LWEParameters, success_probability=0.99):
        """
        Estimate cost of distinguishing a 0-dimensional LWE instance from uniformly random,
        which is essentially the number of samples required.

        :param params: LWE parameters
        :param success_probability: the targeted success probability
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: memory requirement in integers mod q.
        - ``m``: Required number of samples to distinguish.

        EXAMPLE::

            >>> from estimator import *
            >>> from estimator.lwe_guess import distinguish
            >>> params = LWE.Parameters(n=0, q=2 ** 32, Xs=ND.Binary, Xe=ND.DiscreteGaussian(2 ** 32))
            >>> distinguish(params)
            rop: ≈2^60.0, mem: ≈2^60.0, m: ≈2^60.0

        """

        if params.n > 0:
            raise OutOfBoundsError(
                "Secret dimension should be 0 for distinguishing. Try exhaustive search for n > 0."
            )
        m = amplify_sigma(success_probability, sigmaf(params.Xe.stddev), params.q)
        if m > params.m:
            raise InsufficientSamplesError(
                "Not enough samples to distinguish with target advantage."
            )
        return Cost(rop=m, mem=m, m=m).sanity_check()

    __name__ = "distinguish"


distinguish = Distinguisher()
