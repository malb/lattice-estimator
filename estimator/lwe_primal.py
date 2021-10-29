# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE.
"""
from functools import partial

from sage.all import oo, ceil, sqrt, log, RR, ZZ, binomial
from .reduction import BKZ
from .util import binary_search
from .cost import Cost
from .lwe import LWEParameters
from .simulator import Simulator
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai


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
        m = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m)
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
        reduction_cost_model=BKZ.default,
        bkz_model="gsa",
        **kwds,
    ):
        """
        Estimate cost of solving LWE via uSVP reduction.

        :param params: LWE parameters
        :param reduction_cost_model: How to cost BKZ
        :param bkz_model: How to model the shape of a BKZ reduced basis

        EXAMPLE::

            sage: from estimator import *
            sage: print(primal_usvp(Kyber512))
            rop: ≈2^140.9, red: ≈2^140.9, δ: 1.004111, β:  382, d:  973, tag: usvp

            sage: params = LWEParameters(n=200, q=127, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            sage: print(primal_usvp(params, bkz_model="cn11"))
            rop: ≈2^89.0, red: ≈2^89.0, δ: 1.006114, β:  209, d:  388, tag: usvp

            sage: print(primal_usvp(params, bkz_model=Simulator.CN11))
            rop: ≈2^89.0, red: ≈2^89.0, δ: 1.006114, β:  209, d:  388, tag: usvp

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
            f = self.cost_gsa
            start = 40
            stop = 2 * params.n
        elif str(bkz_model).lower() != "gsa":
            try:
                bkz_model = getattr(Simulator, str(bkz_model).upper())
            except AttributeError:
                pass
            cost_gsa = self(
                params,
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
            m=m,
            reduction_cost_model=reduction_cost_model,
            **kwds,
        )

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
    def cost(
        beta: int,
        params: LWEParameters,
        tau: int = 0,
        babai=True,
        mitm=False,
        simulator=Simulator.GSA,
        reduction_cost_model=BKZ.default,
    ):
        """
        Cost of the hybrid attack.

        :param tau: guessing dimension τ
        :param mitm: simulate MITM approach (√ of search space)
        """
        h = len(params.Xs) * params.Xs.density

        d = (params.m + params.n if params.Xs <= params.Xe else params.m) - tau + 1
        scale = PrimalUSVP._scale_factor(params.Xs, params.Xe)

        # 1. Simulate BKZ-β
        r = simulator(d, params.n - tau, params.q, beta, scale=scale)
        bkz_cost = BKZ.cost(reduction_cost_model, beta, d)

        # 2. Required SVP dimension η
        if babai:
            eta = 2
            svp_cost = PrimalHybrid.babai_cost(d)
        else:
            # we scaled the lattice so that χ_e is what we want
            eta = PrimalHybrid.svp_dimension(r, params.Xe)
            svp_cost = BKZ.cost(reduction_cost_model, eta, eta)
            svp_cost["rop"] += PrimalHybrid.babai_cost(d - eta)["rop"]

        # 3. Search
        # We need to do one BDD call at least
        search_space, probability, hw = ZZ(1), 1.0, 0

        # MITM or no MITM
        ssf = sqrt if mitm else lambda x: x

        if tau:
            probability = prob_drop(params.n, h, tau)
            hw = 1
            while hw < h and hw < tau:
                new_search_space = search_space + binomial(tau, hw) * 2 ** hw
                if svp_cost.repeat(ssf(new_search_space))["rop"] < bkz_cost["rop"]:
                    search_space = new_search_space
                    probability += prob_drop(params.n, h, tau, fail=hw)
                    hw += 1
                else:
                    break

            svp_cost = svp_cost.repeat(ssf(search_space))

        if babai is True:
            probability *= prob_babai(r, sqrt(d) * params.Xe.stddev)

        ret = Cost()
        ret["rop"] = bkz_cost["rop"] + svp_cost["rop"]
        ret["pre"] = bkz_cost["rop"]
        ret["svp"] = svp_cost["rop"]
        ret["beta"] = beta
        ret["eta"] = eta
        ret["|S|"] = search_space
        ret["d"] = d
        ret["prob"] = probability

        # 4. Repeat whole experiment ~1/prob times
        ret = ret.repeat(
            prob_amplify(0.99, probability),
            select={
                "rop": True,
                "pre": True,
                "svp": True,
                "beta": False,
                "eta": False,
                "d": False,
                "|S|": False,
                "scale": False,
                "prob": False,
                "pp": False,
            },
        )

        return ret

    def __call__(
        self,
        params: LWEParameters,
        babai: bool = False,
        tau: int = None,
        mitm: bool = True,
        bkz_model="gsa",
        reduction_cost_model=BKZ.default,
        **kwds,
    ):

        try:
            bkz_model = getattr(Simulator, str(bkz_model).upper())
        except AttributeError:
            pass

        if babai is False and tau == 0:
            cost = primal_usvp(
                params,
                bkz_model=bkz_model,
                reduction_cost_model=reduction_cost_model,
                **kwds,
            )

            for b in range(40, cost["beta"])[::-1]:
                cost_curr = self.cost(
                    b,
                    params,
                    tau=tau,
                    babai=babai,
                    mitm=mitm,
                    simulator=bkz_model,
                    reduction_cost_model=reduction_cost_model,
                )
                if cost_curr["rop"] < cost["rop"]:
                    cost = cost_curr
                else:
                    break
            cost["tag"] = cost.data.get("tag", "bdd")
            del cost.data["|S|"]
            del cost.data["prob"]
            del cost.data["repeat"]
            return cost
        else:
            raise NotImplementedError


primal_bdd = partial(PrimalHybrid(), tau=0, mitm=False, babai=False)
