# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

We construct an example LWE instance::

    sage: from estimator import *
    sage: params = LWEParameters(n=384, q=7981, Xs=ND.SparseTernary(384, 16), Xe=ND.CentredBinomial(4))
    sage: params
    LWEParameters(n=384, q=7981, Xs=D(σ=0.29, μ=0.00, n=384), Xe=D(σ=1.41, μ=0.00), m=+Infinity, tag=None)

The simplest (and quickest to estimate) model is solving via uSVP and assuming the Geometric Series
Assumption (GSA)::

    sage: primal_usvp(params, bkz_model="gsa")
    rop: ≈2^86.5, red: ≈2^86.5, δ: 1.006322, β: 198, d: 642, tag: usvp

We get a similar result if we use the ``GSA`` simulator. We do not get the identical result because
we optimize β and d separately::

    sage: primal_usvp(params, bkz_model=Simulator.GSA)
    rop: ≈2^87.3, red: ≈2^87.3, δ: 1.006263, β: 201, d: 603, tag: usvp

To get a more precise answer we may use the CN11 simulator::

    sage: primal_usvp(params, bkz_model=Simulator.CN11)
    rop: ≈2^87.7, red: ≈2^87.7, δ: 1.006244, β: 202, d: 648, tag: usvp

We can then improve on this result by first preprocessing the basis with blocksize β followed by a
single SVP call in dimension η. We call this the BDD approach since this is essentially the same
strategy as preprocessing a basis and then running a CVP solver::

    sage: primal_bdd(params, bkz_model=Simulator.CN11)
    rop: ≈2^83.3, red: ≈2^82.5, svp: ≈2^82.1, β: 184, η: 225, d: 651, tag: bdd

We can improve these results further by exploiting the sparse secret in the hybrid attack, guessing ζ
positions of the secret::

    sage: primal_hybrid(params, bkz_model=Simulator.CN11) # long time
    rop: ≈2^82.8, red: ≈2^82.3, svp: ≈2^81.0, β: 183, η: 186, ζ: 24, |S|: ≈2^20.6, d: 700, prob: 0.991, repeat: 1, ...

"""
from functools import partial
import logging

from sage.all import oo, ceil, sqrt, log, RR, ZZ, binomial, cached_function
from .reduction import BKZ
from .util import binary_search
from .cost import Cost
from .lwe import LWEParameters
from .simulator import Simulator
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai


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
        delta = BKZ.delta(beta)
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
        reduction_cost_model=BKZ.default,
    ):

        delta = BKZ.delta(beta)
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

        return BKZ.cost(reduction_cost_model, beta, d, predicate=lhs <= rhs)

    @staticmethod
    @cached_function
    def cost_simulator(
        beta: int,
        params: LWEParameters,
        simulator,
        m: int = oo,
        tau=None,
        d=None,
        reduction_cost_model=BKZ.default,
    ):
        delta = BKZ.delta(beta)
        if d is None:
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m) + 1
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        tau = params.Xe.stddev if tau is None else tau

        r = simulator(d=d, n=params.n, q=params.q, beta=beta, xi=xi, tau=tau)
        lhs = params.Xe.stddev ** 2 * (beta - 1) + tau ** 2
        if r[d - beta] > lhs:
            cost = BKZ.cost(reduction_cost_model, beta, d)
        else:
            cost = BKZ.cost(reduction_cost_model, beta, d, predicate=False)
        return cost

    def __call__(
        self,
        params: LWEParameters,
        reduction_cost_model=BKZ.default,
        bkz_model="gsa",
        optimize_d=True,
        **kwds,
    ):
        """
        Estimate cost of solving LWE via uSVP reduction.

        :param params: LWE parameters
        :param reduction_cost_model: How to cost BKZ
        :param bkz_model: How to model the shape of a BKZ reduced basis
        :param optimize_d: Attempt to find minimal d, too

        EXAMPLE::

            sage: from estimator import *
            sage: primal_usvp(Kyber512)
            rop: ≈2^140.9, red: ≈2^140.9, δ: 1.004111, β:  382, d:  973, tag: usvp

            sage: params = LWEParameters(n=200, q=127, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            sage: primal_usvp(params, bkz_model="cn11")
            rop: ≈2^89.0, red: ≈2^89.0, δ: 1.006114, β:  209, d:  388, tag: usvp

            sage: primal_usvp(params, bkz_model=Simulator.CN11)
            rop: ≈2^89.0, red: ≈2^89.0, δ: 1.006114, β:  209, d:  388, tag: usvp

            sage: primal_usvp(params, bkz_model=Simulator.CN11, optimize_d=False)
            rop: ≈2^89.1, red: ≈2^89.1, δ: 1.006114, β:  209, d:  400, tag: usvp

        The success condition was formulated in [USENIX:ADPS16]_ and studied/verified in
        [AC:AGVW17,C:DDGR20,PKC:PosVir21]_. The treatment of small secrets is from
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

        if bkz_model == "gsa":
            cost = binary_search(
                self.cost_gsa,
                start=40,
                stop=2 * params.n,
                param="beta",
                predicate=lambda x, best: x["red"] <= best["red"],
                params=params,
                m=m,
                reduction_cost_model=reduction_cost_model,
                **kwds,
            )
            cost["tag"] = "usvp"
            return cost

        try:
            bkz_model = getattr(Simulator, str(bkz_model).upper())
        except AttributeError:
            pass

        # step 0. establish baseline
        cost_gsa = self(
            params,
            reduction_cost_model=reduction_cost_model,
            bkz_model="gsa",
        )
        logging.getLogger("primal").info(f"U0: {repr(cost_gsa)}")

        f = partial(
            self.cost_simulator,
            simulator=bkz_model,
            reduction_cost_model=reduction_cost_model,
            m=m,
            params=params,
        )

        # step 1. find β
        cost = binary_search(
            f,
            param="beta",
            start=cost_gsa["beta"] - 32,
            stop=cost_gsa["beta"] + 128,
            predicate=lambda x, best: x["rop"] <= best["rop"],
            **kwds,
        )
        logging.getLogger("primal").info(f"U1: {repr(cost)}")

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
            logging.getLogger("primal").info(f"U2: {repr(cost)}")

        cost["tag"] = "usvp"
        return cost

    def __repr__(self):
        return "primal_usvp"


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
        simulator=Simulator.GSA,
        reduction_cost_model=BKZ.default,
    ):
        """
        Cost of the hybrid attack.

        :param zeta: guessing dimension ζ
        :param mitm: simulate MITM approach (√ of search space)
        """
        h = ceil(len(params.Xs) * params.Xs.density)

        if d is None:
            delta = BKZ.delta(beta)
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m) + 1
        d -= zeta
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)

        # 1. Simulate BKZ-β
        r = simulator(d, params.n - zeta, params.q, beta, xi=xi)
        bkz_cost = BKZ.cost(reduction_cost_model, beta, d)

        # 2. Required SVP dimension η
        if babai:
            eta = 2
            svp_cost = PrimalHybrid.babai_cost(d)
        else:
            # we xid the lattice so that χ_e is what we want
            eta = PrimalHybrid.svp_dimension(r, params.Xe)
            svp_cost = BKZ.cost(reduction_cost_model, eta, eta)
            svp_cost["rop"] += PrimalHybrid.babai_cost(d - eta)["rop"]

        # 3. Search
        # We need to do one BDD call at least
        search_space, probability, hw = ZZ(1), 1.0, 0

        # MITM or no MITM
        def ssf(x):
            if mitm:
                return RR(sqrt(x))
            else:
                return x

        # e.g. (-1, 1) -> two non-zero per entry
        base = params.Xs.bounds[1] - params.Xs.bounds[0] - 1

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
        Cost.register_impermanent(
            {"|S|": False},
            rop=True,
            red=True,
            svp=True,
            eta=False,
            zeta=False,
            prob=False,
        )
        ret = ret.repeat(
            prob_amplify(0.99, probability),
        )

        return ret

    def cost_zeta(
        self,
        params: LWEParameters,
        bkz_model=Simulator.GSA,
        reduction_cost_model=BKZ.default,
        zeta: int = 0,
        m: int = oo,
        babai: bool = False,
        mitm: bool = True,
        optimize_d=True,
        **kwds,
    ):

        # step 0. establish baseline
        cost = primal_usvp(
            params,
            bkz_model=bkz_model,
            reduction_cost_model=reduction_cost_model,
            optimize_d=False,
            **kwds,
        )
        logging.getLogger("primal").info(f"H0: {repr(cost)}")

        f = partial(
            self.cost,
            params=params,
            zeta=zeta,
            babai=babai,
            mitm=mitm,
            simulator=bkz_model,
            reduction_cost_model=reduction_cost_model,
            m=m,
        )

        # step 1. optimize β
        for b in range(40, cost["beta"])[::-1]:
            cost_curr = f(b)
            if cost_curr["rop"] < cost["rop"]:
                cost = cost_curr
            else:
                break

        logging.getLogger("primal").info(f"H1: {repr(cost)}")

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
            logging.getLogger("primal").info(f"H2: {repr(cost)}")

        if zeta == 0:
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

    def __call__(
        self,
        params: LWEParameters,
        babai: bool = False,
        zeta: int = None,
        mitm: bool = True,
        bkz_model="gsa",
        reduction_cost_model=BKZ.default,
        **kwds,
    ):

        params = LWEParameters.normalize(params)

        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

        try:
            bkz_model = getattr(Simulator, str(bkz_model).upper())
        except AttributeError:
            pass

        f = partial(
            self.cost_zeta,
            params=params,
            bkz_model=bkz_model,
            reduction_cost_model=reduction_cost_model,
            babai=babai,
            mitm=mitm,
            m=m,
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
                return min(cost, cost_0)
            else:
                return f(zeta=zeta)
        else:
            raise NotImplementedError


primal_bdd = partial(PrimalHybrid(), zeta=0, mitm=False, babai=False)
primal_hybrid = partial(PrimalHybrid(), zeta=None, mitm=True, babai=False)
