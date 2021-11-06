# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from functools import partial

from sage.all import oo, ceil, sqrt, log, RR, ZZ, binomial, cached_function
import estimator.reduction as RC
from .util import binary_search
from .cost import Cost
from .lwe import LWEParameters
import estimator.simulator as Simulator
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai
from .io import Logging
from .conf import red_cost_model_default, red_shape_model_default, red_simulator_default


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
        delta = RC.delta(beta)
        a = -log(delta)
        C = log(params.Xe.stddev ** 2 * (beta - 1) + tau ** 2) / 2.0
        b = log(delta) * (2 * beta - 1) + log(params.q) - C
        c = log(tau) + params.n * log(xi) - (params.n + 1) * log(params.q)
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
    ):

        delta = RC.delta(beta)
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        m = min(2 * ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        tau = params.Xe.stddev if tau is None else tau
        d = PrimalUSVP._solve_for_d(params, m, beta, tau, xi) if d is None else d
        assert d <= m + 1

        lhs = log(sqrt(params.Xe.stddev ** 2 * (beta - 1) + tau ** 2))
        rhs = RR(
            log(delta) * (2 * beta - d - 1)
            + (log(tau) + log(xi) * params.n + log(params.q) * (d - params.n - 1)) / d
        )

        return RC.cost(red_cost_model, beta, d, predicate=lhs <= rhs)

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
    ):
        delta = RC.delta(beta)
        if d is None:
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m) + 1
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        tau = params.Xe.stddev if tau is None else tau

        r = simulator(d=d, n=params.n, q=params.q, beta=beta, xi=xi, tau=tau)
        lhs = params.Xe.stddev ** 2 * (beta - 1) + tau ** 2
        if r[d - beta] > lhs:
            cost = RC.cost(red_cost_model, beta, d)
        else:
            cost = RC.cost(red_cost_model, beta, d, predicate=False)
        return cost

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

        :param params: LWE parameters
        :param red_cost_model: How to cost lattice reduction
        :param red_shape_model: How to model the shape of a reduced basis
        :param optimize_d: Attempt to find minimal d, too
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction.
        - ``δ``: Root-Hermite factor targeted by lattice reduction.
        - ``β``: BKZ block size.
        - ``d``: Lattice dimension.

        EXAMPLE::

            >>> from estimator import *
            >>> primal_usvp(Kyber512)
            rop: ≈2^141.2, red: ≈2^141.2, δ: 1.004111, β: 382, d: 973, tag: usvp

            >>> params = LWEParameters(n=200, q=127, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            >>> primal_usvp(params, red_shape_model="cn11")
            rop: ≈2^91.2, red: ≈2^91.2, δ: 1.006114, β: 209, d: 388, tag: usvp

            >>> primal_usvp(params, red_shape_model=Simulator.CN11)
            rop: ≈2^91.2, red: ≈2^91.2, δ: 1.006114, β: 209, d: 388, tag: usvp

            >>> primal_usvp(params, red_shape_model=Simulator.CN11, optimize_d=False)
            rop: ≈2^91.3, red: ≈2^91.3, δ: 1.006114, β: 209, d: 400, tag: usvp

        The success condition was formulated in [USENIX:ADPS16]_ and studied/verified in
        [AC:AGVW17]_, [C:DDGR20]_, [PKC:PosVir21]_. The treatment of small secrets is from
        [ACISP:BaiGal14]_.

        .. [ACISP:BaiGal14] Bai, S., & Galbraith, S. D. (2014). Lattice decoding attacks on binary
           LWE. In W. Susilo, & Y. Mu, ACISP 14 (pp. 322–337). : Springer, Heidelberg.

        .. [USENIX:ADPS16] Alkim, E., L\'eo Ducas, Thomas P\"oppelmann, & Schwabe, P. (2016).
           Post-quantum key exchange - A new hope. In T. Holz, & S. Savage, USENIX Security 2016 (pp.
           327–343). : USENIX Association.

        .. [AC:AGVW17] Albrecht, M. R., Florian Göpfert, Virdia, F., & Wunderer, T. (2017).
           Revisiting the expected cost of solving uSVP and applications to LWE. In T. Takagi, & T.
           Peyrin, ASIACRYPT 2017, Part I (pp. 297–322). : Springer, Heidelberg.

        .. [C:DDGR20] Dana Dachman-Soled, Léo Ducas, Gong, H., & M\'elissa Rossi (2020). LWE with
           side information: Attacks and concrete security estimation. In D. Micciancio, & T.
           Ristenpart, CRYPTO~2020, Part~II (pp. 329–358). : Springer, Heidelberg.

        .. [PKC:PosVir21] Postlethwaite, E. W., & Virdia, F. (2021). On the success probability of
           solving unique SVP via BKZ. In J. Garay, PKC 2021, Part I (pp. 68–98). : Springer,
           Heidelberg.
        """
        params = LWEParameters.normalize(params)

        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

        if red_shape_model == "gsa":
            cost = binary_search(
                self.cost_gsa,
                start=40,
                stop=2 * params.n,
                param="beta",
                predicate=lambda x, best: x["rop"] <= best["rop"],
                params=params,
                m=m,
                red_cost_model=red_cost_model,
                **kwds,
            )
            cost["tag"] = "usvp"
            return cost

        try:
            red_shape_model = getattr(Simulator, str(red_shape_model).upper())
        except AttributeError:
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
        cost = binary_search(
            f,
            param="beta",
            start=cost_gsa["beta"] - ceil(0.10 * cost_gsa["beta"]),
            stop=cost_gsa["beta"] + ceil(0.20 * cost_gsa["beta"]),
            predicate=lambda x, best: x["rop"] <= best["rop"],
            **kwds,
        )
        Logging.log("usvp", log_level, f"Opt-β: {repr(cost)}")

        if cost and optimize_d:
            # step 2. find d
            cost = binary_search(
                f,
                param="d",
                start=params.n,
                stop=cost["d"],
                predicate=lambda x, best: x["rop"] <= best["rop"],
                beta=cost["beta"],
                **kwds,
            )
            Logging.log("usvp", log_level + 1, f"Opt-d: {repr(cost)}")

        cost["tag"] = "usvp"
        return cost

    __name__ = "primal_usvp"


primal_usvp = PrimalUSVP()


class PrimalHybrid:
    @classmethod
    def babai_cost(cls, d):
        return Cost(rop=d ** 2)

    @classmethod
    def svp_dimension(cls, r, D):
        """
        Return η for a given lattice shape and distance.

        :param r: squared Gram-Schmidt norms

        """
        from fpylll.util import gaussian_heuristic

        d = len(r)
        for i, _ in enumerate(r):
            if gaussian_heuristic(r[i:]) < D.stddev ** 2 * (d - i):
                return ZZ(d - (i - 1))
        return ZZ(2)

    @staticmethod
    @cached_function
    def cost(
        beta: int,
        params: LWEParameters,
        zeta: int = 0,
        babai=True,
        mitm=False,
        m: int = oo,
        d: int = None,
        simulator=red_simulator_default,
        red_cost_model=red_cost_model_default,
        log_level=5,
    ):
        """
        Cost of the hybrid attack.

        :param zeta: guessing dimension ζ ≥ 0
        :param babai: insist on Babai's algorithm for finding close vectors
        :param mitm: simulate MITM approach (√ of search space)
        :param m: We accept the number of samples to consider from the calling function
        :param d: We optionally accept the dimension to pick

        .. note :: This is the lowest level function that runs no optimization, it merely reports
           costs.

        """
        # the number of non-zero entries
        h = ceil(len(params.Xs) * params.Xs.density)

        if d is None:
            delta = RC.delta(beta)
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m) + 1
        d -= zeta
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)

        # 1. Simulate BKZ-β
        # TODO: pick τ
        r = simulator(d, params.n - zeta, params.q, beta, xi=xi)
        bkz_cost = RC.cost(red_cost_model, beta, d)

        # 2. Required SVP dimension η
        if babai:
            eta = 2
            svp_cost = PrimalHybrid.babai_cost(d)
        else:
            # we scaled the lattice so that χ_e is what we want
            eta = PrimalHybrid.svp_dimension(r, params.Xe)
            svp_cost = RC.cost(red_cost_model, eta, eta)
            # when η ≪ β, lifting may be a bigger cost
            svp_cost["rop"] += PrimalHybrid.babai_cost(d - eta)["rop"]

        # 3. Search
        # We need to do one BDD call at least
        search_space, probability, hw = ZZ(1), 1.0, 0

        # MITM or no MITM
        def ssf(x):
            # TODO: this is rather clumsy as a model
            if mitm:
                return RR(sqrt(x))
            else:
                return x

        # e.g. (-1, 1) -> two non-zero per entry
        base = params.Xs.bounds[1] - params.Xs.bounds[0]

        if zeta:
            probability = RR(prob_drop(params.n, h, zeta))
            hw = 1
            while hw < h and hw < zeta and base < oo:
                new_search_space = search_space + binomial(zeta, hw) * base ** hw
                if svp_cost.repeat(ssf(new_search_space))["rop"] < bkz_cost["rop"]:
                    search_space = new_search_space
                    probability += prob_drop(params.n, h, zeta, fail=hw)
                    hw += 1
                else:
                    break

            svp_cost = svp_cost.repeat(ssf(search_space))

        if babai is True:
            probability *= RR(prob_babai(r, sqrt(d) * params.Xe.stddev))

        Cost.register_impermanent(
            {"|S|": False},
            rop=True,
            red=True,
            svp=True,
            eta=False,
            zeta=False,
            prob=False,
        )

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

        # 4. Repeat whole experiment ~1/prob times
        ret = ret.repeat(
            prob_amplify(0.99, probability),
        )

        return ret

    def cost_zeta(
        self,
        zeta: int,
        params: LWEParameters,
        red_shape_model=red_simulator_default,
        red_cost_model=red_cost_model_default,
        m: int = oo,
        babai: bool = False,
        mitm: bool = True,
        optimize_d=True,
        baseline_cost=None,
        log_level=5,
        **kwds,
    ):
        """
        This function optimizes costs for a fixed guessing dimension ζ.
        """

        # step 0. establish baseline
        if baseline_cost is None:
            baseline_cost = primal_usvp(
                params,
                red_shape_model=red_shape_model,
                red_cost_model=red_cost_model,
                optimize_d=False,
                log_level=log_level + 1,
                **kwds,
            )
            Logging.log("bdd", log_level + 1, f"H0: {repr(baseline_cost)}")

        f = partial(
            self.cost,
            params=params,
            zeta=zeta,
            babai=babai,
            mitm=mitm,
            simulator=red_shape_model,
            red_cost_model=red_cost_model,
            m=m,
        )

        # step 1. optimize β
        for b in range(40, baseline_cost["beta"])[::-1]:
            cost = f(b)
            if cost["rop"] < baseline_cost["rop"]:
                baseline_cost = cost
            else:
                break

        Logging.log("bdd", log_level, f"H1: {repr(cost)}")

        # step 2. optimize d
        if cost and cost.get("tag", "XXX") != "usvp" and optimize_d:
            cost = binary_search(
                f,
                param="d",
                start=params.n,
                stop=cost["d"],
                predicate=lambda x, best: x["rop"] <= best["rop"],
                beta=cost["beta"],
                **kwds,
            )
            Logging.log("bdd", log_level + 1, f"H2: {repr(cost)}")

        return cost

    def __call__(
        self,
        params: LWEParameters,
        babai: bool = False,
        zeta: int = None,
        mitm: bool = True,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=1,
        **kwds,
    ):

        params = LWEParameters.normalize(params)

        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

        try:
            red_shape_model = getattr(Simulator, str(red_shape_model).upper())
        except AttributeError:
            pass

        f = partial(
            self.cost_zeta,
            params=params,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            babai=babai,
            mitm=mitm,
            m=m,
            log_level=log_level,
            **kwds,
        )

        if babai is False:
            if zeta is None:
                cost_0 = f(zeta=0)
                cost = binary_search(
                    f,
                    start=0,
                    stop=params.n,
                    param="zeta",
                    predicate=lambda x, best: x <= best,
                    optimize_d=False,
                )
                cost = min(cost, cost_0)
            else:
                cost = f(zeta=zeta)
        else:
            raise NotImplementedError

        if cost.get("zeta", 0) == 0:
            cost["tag"] = cost.get("tag", "bdd")
            try:
                del cost["|S|"]
                del cost["prob"]
                del cost["repeat"]
                del cost["zeta"]
            except KeyError:
                pass
        else:
            cost["tag"] = cost.get("tag", "hybrid")
        return cost

    __name__ = "primal_hybrid"


primal_hybrid = PrimalHybrid()


def primal_bdd(
    params: LWEParameters,
    red_shape_model=red_shape_model_default,
    red_cost_model=red_cost_model_default,
    log_level=1,
    **kwds,
):
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
