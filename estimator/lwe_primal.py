# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from functools import partial

from sage.all import oo, ceil, sqrt, log, RR, ZZ, binomial, cached_function
from .reduction import delta as deltaf
from .reduction import cost as costf
from .util import local_minimum
from .cost import Cost
from .lwe_parameters import LWEParameters
from .simulator import normalize as simulator_normalize
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai
from .prob import mitm_babai_probability
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .conf import red_simulator as red_simulator_default


class PrimalUSVP:
    """
    Estimate cost of solving LWE via uSVP reduction.
    """

    @staticmethod
    def _xi_factor(Xs, Xe):
        xi = RR(1)
        if Xs < Xe:
            xi = Xe.stddev / Xs.stddev
        return xi

    @staticmethod
    def _solve_for_d(params, m, beta, tau, xi):
        """
        Find smallest d ∈ [n,m] to satisfy uSVP condition.

        If no such d exists, return the upper bound m.
        """
        # Find the smallest d ∈ [n,m] s.t. a*d^2 + b*d + c >= 0
        delta = deltaf(beta)
        a = -log(delta)

        if not tau:
            C = log(params.Xe.stddev**2 * (beta - 1)) / 2.0
            c = params.n * log(xi) - (params.n + 1) * log(params.q)

        else:
            C = log(params.Xe.stddev**2 * (beta - 1) + tau**2) / 2.0
            c = log(tau) + params.n * log(xi) - (params.n + 1) * log(params.q)

        b = log(delta) * (2 * beta - 1) + log(params.q) - C
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
    @cached_function
    def cost_gsa(
        beta: int,
        params: LWEParameters,
        m: int = oo,
        tau=None,
        d=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
    ):
        delta = deltaf(beta)
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        m = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        tau = params.Xe.stddev if tau is None else tau
        # Account for homogeneous instances
        if params._homogeneous:
            tau = False  # Tau false ==> instance is homogeneous

        d = PrimalUSVP._solve_for_d(params, m, beta, tau, xi) if d is None else d
        if d < beta:
            d = beta
        # if d == β we assume one SVP call, otherwise poly calls. This makes the cost curve jump, so
        # we avoid it here.
        if d == beta and d < m:
            d += 1
        assert d <= m + 1

        if not tau:
            lhs = log(sqrt(params.Xe.stddev**2 * (beta - 1)))
            rhs = RR(
                log(delta) * (2 * beta - d - 1)
                + (log(xi) * params.n + log(params.q) * (d - params.n - 1)) / d
            )

        else:
            lhs = log(sqrt(params.Xe.stddev**2 * (beta - 1) + tau**2))
            rhs = RR(
                log(delta) * (2 * beta - d - 1)
                + (log(tau) + log(xi) * params.n + log(params.q) * (d - params.n - 1)) / d
            )

        return costf(red_cost_model, beta, d, predicate=lhs <= rhs)

    @staticmethod
    @cached_function
    def cost_simulator(
        beta: int,
        params: LWEParameters,
        simulator,
        m: int = oo,
        tau=None,
        d=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
    ):
        delta = deltaf(beta)
        if d is None:
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m) + 1
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        tau = params.Xe.stddev if tau is None else tau

        if params._homogeneous:
            tau = False
            d -= 1  # Remove extra dimension in homogeneous instances

        r = simulator(d=d, n=params.n, q=params.q, beta=beta, xi=xi, tau=tau)

        if not tau:
            lhs = params.Xe.stddev**2 * (beta - 1)

        else:
            lhs = params.Xe.stddev**2 * (beta - 1) + tau**2

        predicate = r[d - beta] > lhs

        return costf(red_cost_model, beta, d, predicate=predicate)

    def __call__(
        self,
        params: LWEParameters,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        optimize_d=True,
        log_level=1,
        **kwds,
    ):
        """
        Estimate cost of solving LWE via uSVP reduction.

        :param params: LWE parameters.
        :param red_cost_model: How to cost lattice reduction.
        :param red_shape_model: How to model the shape of a reduced basis.
        :param optimize_d: Attempt to find minimal d, too.
        :return: A cost dictionary.

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction.
        - ``δ``: Root-Hermite factor targeted by lattice reduction.
        - ``β``: BKZ block size.
        - ``d``: Lattice dimension.

        EXAMPLE::

            >>> from estimator import *
            >>> LWE.primal_usvp(schemes.Kyber512)
            rop: ≈2^143.8, red: ≈2^143.8, δ: 1.003941, β: 406, d: 998, tag: usvp

            >>> params = LWE.Parameters(n=200, q=127, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            >>> LWE.primal_usvp(params, red_shape_model="cn11")
            rop: ≈2^87.6, red: ≈2^87.6, δ: 1.006114, β: 209, d: 388, tag: usvp

            >>> LWE.primal_usvp(params, red_shape_model=Simulator.CN11)
            rop: ≈2^87.6, red: ≈2^87.6, δ: 1.006114, β: 209, d: 388, tag: usvp

            >>> LWE.primal_usvp(params, red_shape_model=Simulator.CN11, optimize_d=False)
            rop: ≈2^87.6, red: ≈2^87.6, δ: 1.006114, β: 209, d: 400, tag: usvp

            >>> params = LWE.Parameters(n=384, q=2**7, Xs=ND.Uniform(0, 1), Xe=ND.CenteredBinomial(8), m=2*384)
            >>> LWE.primal_usvp(params, red_cost_model=RC.BDGL16)  # Issue #87
            rop: ≈2^161.8, red: ≈2^161.8, δ: 1.003634, β: 456, d: 595, tag: usvp

            >>> Xe=ND.DiscreteGaussian(stddev=3.19)
            >>> params = LWE.Parameters(n=1030, m=2060, q=2**64, Xs=ND.Uniform(0, 1), Xe=Xe)
            >>> LWE.primal_usvp(params, red_cost_model=RC.BDGL16)  # Issue 95
            rop: ≈2^56.6, red: ≈2^56.6, δ: 1.009686, β: 91, d: 1618, tag: usvp

        The success condition was formulated in [USENIX:ADPS16]_ and studied/verified in
        [AC:AGVW17]_, [C:DDGR20]_, [PKC:PosVir21]_. The treatment of small secrets is from
        [ACISP:BaiGal14]_.

        """
        params = LWEParameters.normalize(params)
        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m
        if red_shape_model == "gsa":
            with local_minimum(40, max(2 * params.n, 41), precision=5) as it:
                for beta in it:
                    cost = self.cost_gsa(
                        beta=beta, params=params, m=m, red_cost_model=red_cost_model, **kwds
                    )
                    it.update(cost)
                for beta in it.neighborhood:
                    cost = self.cost_gsa(
                        beta=beta, params=params, m=m, red_cost_model=red_cost_model, **kwds
                    )
                    it.update(cost)
                cost = it.y
            cost["tag"] = "usvp"
            cost["problem"] = params
            return cost.sanity_check()

        try:
            red_shape_model = simulator_normalize(red_shape_model)
        except ValueError:
            pass

        # step 0. establish baseline
        cost_gsa = self(
            params,
            red_cost_model=red_cost_model,
            red_shape_model="gsa",
        )

        Logging.log("usvp", log_level + 1, f"GSA: {repr(cost_gsa)}")

        f = partial(
            self.cost_simulator,
            simulator=red_shape_model,
            red_cost_model=red_cost_model,
            m=m,
            params=params,
        )

        # step 1. find β

        with local_minimum(
            max(cost_gsa["beta"] - ceil(0.10 * cost_gsa["beta"]), 40),
            max(cost_gsa["beta"] + ceil(0.20 * cost_gsa["beta"]), 40),
        ) as it:
            for beta in it:
                it.update(f(beta=beta, **kwds))
            cost = it.y

        Logging.log("usvp", log_level, f"Opt-β: {repr(cost)}")

        if cost and optimize_d:
            # step 2. find d
            with local_minimum(params.n, stop=cost["d"] + 1) as it:
                for d in it:
                    it.update(f(d=d, beta=cost["beta"], **kwds))
                cost = it.y
            Logging.log("usvp", log_level + 1, f"Opt-d: {repr(cost)}")

        cost["tag"] = "usvp"
        cost["problem"] = params
        return cost.sanity_check()

    __name__ = "primal_usvp"


primal_usvp = PrimalUSVP()


class PrimalHybrid:
    @classmethod
    def babai_cost(cls, d):
        return Cost(rop=max(d, 1) ** 2)

    @classmethod
    def svp_dimension(cls, r, D):
        """
        Return η for a given lattice shape and distance.

        :param r: squared Gram-Schmidt norms

        """
        from math import lgamma, log, exp, pi

        def ball_log_vol(n):
            return (n / 2.0) * log(pi) - lgamma(n / 2.0 + 1)

        def gaussian_heuristic_log_input(r):
            n = len(list(r))
            log_vol = sum(r)
            log_gh = 1.0 / n * (log_vol - 2 * ball_log_vol(n))
            return exp(log_gh)

        d = len(r)
        r = [log(x) for x in r]

        if d > 4096:
            for i, _ in enumerate(r):
                # chosen since RC.ADPS16(1754, 1754).log(2.) = 512.168000000000
                j = d - 1754 + i
                if (j < d) and (gaussian_heuristic_log_input(r[j:]) < D.stddev**2 * (d - j)):
                    return ZZ(d - (j - 1))
            return ZZ(2)

        else:
            for i, _ in enumerate(r):
                if gaussian_heuristic_log_input(r[i:]) < D.stddev**2 * (d - i):
                    return ZZ(d - (i - 1))
            return ZZ(2)

    @staticmethod
    @cached_function
    def cost(
        beta: int,
        params: LWEParameters,
        zeta: int = 0,
        babai=False,
        mitm=False,
        m: int = oo,
        d: int = None,
        simulator=red_simulator_default,
        red_cost_model=red_cost_model_default,
        log_level=5,
    ):
        """
        Cost of the hybrid attack.

        :param beta: Block size.
        :param params: LWE parameters.
        :param zeta: Guessing dimension ζ ≥ 0.
        :param babai: Insist on Babai's algorithm for finding close vectors.
        :param mitm: Simulate MITM approach (√ of search space).
        :param m: We accept the number of samples to consider from the calling function.
        :param d: We optionally accept the dimension to pick.

        .. note :: This is the lowest level function that runs no optimization, it merely reports
           costs.

        """
        if d is None:
            delta = deltaf(beta)
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m) + 1
        d -= zeta

        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        tau = 1
        # 1. Simulate BKZ-β
        # TODO: pick τ as non default value

        if params._homogeneous:
            tau = False
            d -= 1

        r = simulator(d, params.n - zeta, params.q, beta, xi=xi, tau=tau, dual=True)

        bkz_cost = costf(red_cost_model, beta, d)

        # 2. Required SVP dimension η
        if babai:
            eta = 2
            svp_cost = PrimalHybrid.babai_cost(d)
        else:
            # we scaled the lattice so that χ_e is what we want
            eta = PrimalHybrid.svp_dimension(r, params.Xe)
            if eta > d:
                # Lattice reduction was not strong enough to "reveal" the LWE solution.
                # A larger `beta` should perhaps be attempted.
                return Cost(rop=oo)
            svp_cost = costf(red_cost_model, eta, eta)
            # when η ≪ β, lifting may be a bigger cost
            svp_cost["rop"] += PrimalHybrid.babai_cost(d - eta)["rop"]

        # 3. Search
        # We need to do one BDD call at least
        search_space, probability, hw = 1, 1.0, 0

        # MITM or no MITM
        # TODO: this is rather clumsy as a model
        def ssf(x):
            if mitm:
                return RR(sqrt(x))
            else:
                return x

        # e.g. (-1, 1) -> two non-zero per entry
        base = params.Xs.bounds[1] - params.Xs.bounds[0]

        if zeta:
            # the number of non-zero entries
            h = params.Xs.hamming_weight
            probability = RR(prob_drop(params.n, h, zeta))
            hw = 1
            while hw < min(h, zeta):
                new_search_space = binomial(zeta, hw) * base**hw
                if svp_cost.repeat(ssf(search_space + new_search_space))["rop"] >= bkz_cost["rop"]:
                    break
                search_space += new_search_space
                probability += prob_drop(params.n, h, zeta, fail=hw)
                hw += 1

            svp_cost = svp_cost.repeat(ssf(search_space))

        if mitm and zeta > 0:
            if babai:
                probability *= mitm_babai_probability(r, params.Xe.stddev)
            else:
                # TODO: the probability in this case needs to be analysed
                probability *= 1

        if eta <= 20 and d >= 0:  # NOTE: η: somewhat arbitrary bound, d: we may guess it all
            probability *= RR(prob_babai(r, sqrt(d) * params.Xe.stddev))

        ret = Cost()
        ret["rop"] = bkz_cost["rop"] + svp_cost["rop"]
        ret["red"] = bkz_cost["rop"]
        ret["svp"] = svp_cost["rop"]
        ret["beta"] = beta
        ret["eta"] = eta
        ret["zeta"] = zeta
        ret["|S|"] = search_space
        ret["d"] = d
        ret["prob"] = probability

        ret.register_impermanent(
            {"|S|": False},
            rop=True,
            red=True,
            svp=True,
            eta=False,
            zeta=False,
            prob=False,
        )

        # 4. Repeat whole experiment ~1/prob times
        if probability and not RR(probability).is_NaN():
            ret = ret.repeat(
                prob_amplify(0.99, probability),
            )
        else:
            return Cost(rop=oo)

        return ret

    @classmethod
    def cost_zeta(
        cls,
        zeta: int,
        params: LWEParameters,
        red_shape_model=red_simulator_default,
        red_cost_model=red_cost_model_default,
        m: int = oo,
        babai: bool = True,
        mitm: bool = True,
        optimize_d=True,
        log_level=5,
        **kwds,
    ):
        """
        This function optimizes costs for a fixed guessing dimension ζ.
        """

        # step 0. establish baseline
        baseline_cost = primal_usvp(
            params,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            optimize_d=False,
            log_level=log_level + 1,
            **kwds,
        )
        Logging.log("bdd", log_level, f"H0: {repr(baseline_cost)}")

        f = partial(
            cls.cost,
            params=params,
            zeta=zeta,
            babai=babai,
            mitm=mitm,
            simulator=red_shape_model,
            red_cost_model=red_cost_model,
            m=m,
            **kwds,
        )

        # step 1. optimize β
        with local_minimum(
            40, baseline_cost["beta"] + 1, precision=2, log_level=log_level + 1
        ) as it:
            for beta in it:
                it.update(f(beta))
            for beta in it.neighborhood:
                it.update(f(beta))
            cost = it.y

        Logging.log("bdd", log_level, f"H1: {cost!r}")

        # step 2. optimize d
        if cost and cost.get("tag", "XXX") != "usvp" and optimize_d:
            with local_minimum(
                params.n, cost["d"] + cost["zeta"] + 1, log_level=log_level + 1
            ) as it:
                for d in it:
                    it.update(f(beta=cost["beta"], d=d))
                cost = it.y
            Logging.log("bdd", log_level, f"H2: {cost!r}")

        if cost is None:
            return Cost(rop=oo)
        return cost

    def __call__(
        self,
        params: LWEParameters,
        babai: bool = True,
        zeta: int = None,
        mitm: bool = True,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=1,
        **kwds,
    ):
        """
        Estimate the cost of the hybrid attack and its variants.

        :param params: LWE parameters.
        :param zeta: Guessing dimension ζ ≥ 0.
        :param babai: Insist on Babai's algorithm for finding close vectors.
        :param mitm: Simulate MITM approach (√ of search space).
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction.
        - ``δ``: Root-Hermite factor targeted by lattice reduction.
        - ``β``: BKZ block size.
        - ``η``: Dimension of the final BDD call.
        - ``ζ``: Number of guessed coordinates.
        - ``|S|``: Guessing search space.
        - ``prob``: Probability of success in guessing.
        - ``repeat``: How often to repeat the attack.
        - ``d``: Lattice dimension.

        - When ζ = 0 this function essentially estimates the BDD strategy as given in [RSA:LiuNgu13]_.
        - When ζ ≠ 0 and ``babai=True`` this function estimates the hybrid attack as given in
          [C:HowgraveGraham07]_
        - When ζ ≠ 0 and ``babai=False`` this function estimates the hybrid attack as given in
          [SAC:AlbCurWun19]_

        EXAMPLES::

            >>> from estimator import *
            >>> params = schemes.Kyber512.updated(Xs=ND.SparseTernary(16))
            >>> LWE.primal_hybrid(params, mitm=False, babai=False)
            rop: ≈2^91.5, red: ≈2^90.7, svp: ≈2^90.2, β: 178, η: 21, ζ: 256, |S|: ≈2^56.6, d: 531, prob: 0.003, ↻: 1...

            >>> LWE.primal_hybrid(params, mitm=False, babai=True)
            rop: ≈2^88.7, red: ≈2^88.0, svp: ≈2^87.2, β: 98, η: 2, ζ: 323, |S|: ≈2^39.7, d: 346, prob: ≈2^-28.4, ↻: ...

            >>> LWE.primal_hybrid(params, mitm=True, babai=False)
            rop: ≈2^74.1, red: ≈2^73.7, svp: ≈2^71.9, β: 104, η: 16, ζ: 320, |S|: ≈2^77.1, d: 359, prob: ≈2^-12.3, ↻...

            >>> LWE.primal_hybrid(params, mitm=True, babai=True)
            rop: ≈2^85.8, red: ≈2^84.8, svp: ≈2^84.8, β: 105, η: 2, ζ: 366, |S|: ≈2^85.1, d: 315, prob: ≈2^-23.4, ↻:...

        TESTS:

        We test a trivial instance::

            >>> params = LWE.Parameters(2**10, 2**100, ND.DiscreteGaussian(3.19), ND.DiscreteGaussian(3.19))
            >>> LWE.primal_bdd(params)
            rop: ≈2^43.7, red: ≈2^43.7, svp: ≈2^22.1, β: 40, η: 2, d: 1516, tag: bdd

        We also test a LWE instance with a large error (coming from issue #106)::

            >>> LWE.primal_bdd(LWE.Parameters(n=256, q=12289, Xs=ND.UniformMod(2), Xe=ND.UniformMod(1024)))
            rop: ≈2^116.3, red: ≈2^41.4, svp: ≈2^116.3, β: 40, η: 340, d: 340, tag: bdd

            >>> LWE.primal_bdd(LWE.Parameters(n=700, q=2**64, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2**59)))
            rop: ≈2^263.2, red: ≈2^42.8, svp: ≈2^263.2, β: 40, η: 867, d: 867, tag: bdd


        """

        if zeta == 0:
            tag = "bdd"
        else:
            tag = "hybrid"

        params = LWEParameters.normalize(params)

        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

        red_shape_model = simulator_normalize(red_shape_model)

        f = partial(
            self.cost_zeta,
            params=params,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            babai=babai,
            mitm=mitm,
            m=m,
            log_level=log_level + 1,
        )

        if zeta is None:
            # Find the smallest value for zeta such that the square root of the search space for
            # zeta is larger than the number of operations to solve uSVP on the whole LWE instance
            # (without guessing).
            usvp_cost = primal_usvp(params, red_cost_model=red_cost_model)["rop"]
            zeta_max = params.n
            while zeta_max < params.n and sqrt(params.Xs.resize(zeta_max).support_size()) < usvp_cost:
                zeta_max += 1

            with local_minimum(0, min(zeta_max, params.n), log_level=log_level) as it:
                for zeta in it:
                    it.update(f(zeta=zeta, optimize_d=False, **kwds))
            # TODO: this should not be required
            cost = min(it.y, f(0, optimize_d=False, **kwds))
        else:
            cost = f(zeta=zeta)

        cost["tag"] = tag
        cost["problem"] = params

        if tag == "bdd":
            for k in ("|S|", "prob", "repetitions", "zeta"):
                try:
                    del cost[k]
                except KeyError:
                    pass

        return cost.sanity_check()

    __name__ = "primal_hybrid"


primal_hybrid = PrimalHybrid()


def primal_bdd(
    params: LWEParameters,
    red_shape_model=red_shape_model_default,
    red_cost_model=red_cost_model_default,
    log_level=1,
    **kwds,
):
    """
    Estimate the cost of the BDD approach as given in [RSA:LiuNgu13]_.

    :param params: LWE parameters.
    :param red_cost_model: How to cost lattice reduction
    :param red_shape_model: How to model the shape of a reduced basis

    """

    return primal_hybrid(
        params,
        zeta=0,
        mitm=False,
        babai=False,
        red_shape_model=red_shape_model,
        red_cost_model=red_cost_model,
        log_level=log_level,
        **kwds,
    )
