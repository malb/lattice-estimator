# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from functools import partial

from sage.all import oo, ceil, sqrt, log, RR, ZZ, cached_function
from .reduction import delta as deltaf
from .reduction import cost as costf
from .util import local_minimum
from .cost import Cost
from .lwe_parameters import LWEParameters
from .simulator import normalize as simulator_normalize
from .prob import guessing_set_and_hit_probability
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai
from .prob import mitm_babai_probability
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .conf import max_beta as max_beta_global
from scipy.optimize import minimize_scalar


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

        If no such d exists, return oo.
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
        if disc < 0:  # no solution, return oo
            return oo

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
        return oo

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
        tau = params.Xe.stddev if tau is None else tau
        # Account for homogeneous instances
        if params._homogeneous:
            tau = False  # Tau false ==> instance is homogeneous

        d = PrimalUSVP._solve_for_d(params, m, beta, tau, xi) if d is None else d
        # this beta is not sufficient to reveal the solution with the number of samples m
        if d == oo:
            return Cost(rop=oo)
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
            d = max(d, beta)
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
            rop: ≈2^87.5, red: ≈2^87.5, δ: 1.006114, β: 209, d: 388, tag: usvp

            >>> LWE.primal_usvp(params, red_shape_model=Simulator.CN11)
            rop: ≈2^87.5, red: ≈2^87.5, δ: 1.006114, β: 209, d: 388, tag: usvp

            >>> LWE.primal_usvp(params, red_shape_model=Simulator.CN11, optimize_d=False)
            rop: ≈2^87.6, red: ≈2^87.6, δ: 1.006114, β: 209, d: 400, tag: usvp

            >>> params = LWE.Parameters(n=384, q=2**7, Xs=ND.Uniform(0, 1), Xe=ND.CenteredBinomial(8), m=2*384)
            >>> LWE.primal_usvp(params, red_cost_model=RC.BDGL16)  # Issue #87
            rop: ≈2^161.8, red: ≈2^161.8, δ: 1.003634, β: 456, d: 595, tag: usvp

            >>> Xe=ND.DiscreteGaussian(stddev=3.19)
            >>> params = LWE.Parameters(n=1030, m=2060, q=2**64, Xs=ND.Uniform(0, 1), Xe=Xe)
            >>> LWE.primal_usvp(params, red_cost_model=RC.BDGL16)  # Issue #95
            rop: ≈2^53.1, red: ≈2^53.1, δ: 1.010374, β: 78, d: 1933, tag: usvp

            # small n examples (Issue #181)
            >>> params = LWE.Parameters(n=11, q = 2**128, Xs=ND.UniformMod(2**128), Xe=ND.UniformMod(2**124))
            >>> LWE.primal_usvp(params)
            rop: ≈2^126.0, red: ≈2^126.0, δ: 1.004356, β: 351, d: 455, tag: usvp
            >>> LWE.primal_usvp(params, red_shape_model=Simulator.CN11)
            rop: ≈2^127.1, red: ≈2^127.1, δ: 1.004315, β: 356, d: 443, tag: usvp

        The success condition was formulated in [USENIX:ADPS16]_ and studied/verified in
        [AC:AGVW17]_, [C:DDGR20]_, [PKC:PosVir21]_. The treatment of small secrets is from
        [ACISP:BaiGal14]_.

        """
        params = LWEParameters.normalize(params)

        if params.Xs <= params.Xe:
            # allow for a larger embedding lattice dimension: Bai and Galbraith
            m = params.m + params.n
        else:
            m = params.m
        if red_shape_model == "gsa":
            precision = 5
            max_beta = max(min(max_beta_global, m), 40 + precision)
            with local_minimum(40, max_beta, precision=precision) as it:
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
            with local_minimum(max(params.n, cost["beta"]), stop=cost["d"] + 1) as it:
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
    def svp_dimension(cls, r, D, is_homogeneous=False):
        """
        Return required svp dimension for a given lattice shape and distance.

        :param r: squared Gram-Schmidt norms

        """
        from math import lgamma, log, pi

        def ball_log_vol(n):
            return (n / 2.0) * log(pi) - lgamma(n / 2.0 + 1)

        # If B is a basis with GSO profiles r, this returns an estimate for the shortest vector in the lattice
        # [ B | * ]
        # [ 0 |tau]
        # if the tau is None, the instance is homogeneous, and we omit the final row/column.
        def svp_gaussian_heuristic_log_input(r, tau):
            if tau is None:
                n = len(list(r))
                log_vol = sum(r)
            else:
                n = len(list(r)) + 1
                log_vol = sum(r) + 2 * log(tau)
            log_gh = 1.0 / n * (log_vol - 2 * ball_log_vol(n))
            return log_gh

        d = len(r)
        r = [log(x) for x in r]

        if d > 4096:
            # chosen since RC.ADPS16(1754, 1754).log(2.) = 512.168000000000
            min_i = d - 1754
        else:
            min_i = 0

        if is_homogeneous:
            tau = None
            for i in range(min_i, d):
                if svp_gaussian_heuristic_log_input(r[i:], tau) < log(D.stddev**2 * (d - i)):
                    return ZZ(d - (i - 1))
            return ZZ(2)

        else:
            # we look for the largest i such that (pi_i(e), tau) is shortest in the embedding lattice
            # [pi_i(B) | * ]
            # [   0    |tau]
            tau = D.stddev
            for i in range(min_i, d):
                if svp_gaussian_heuristic_log_input(r[i:], tau) < log(D.stddev**2 * (d - i) + tau ** 2):
                    return ZZ(d - (i - 1) + 1)
            return ZZ(2)

    @classmethod
    def svp_dimension_gsa(cls, d, log_total_vol, log_delta, D, is_homogeneous=False):
        """
        Return required svp dimension assuming the GSA on a lattice with a given volume and rank.

        """
        from math import lgamma, log, pi

        def log_projected_vol(i):
            return (d - i) / d * log_total_vol - i * (d - i) * log_delta

        def ball_log_vol(n):
            return (n / 2.0) * log(pi) - lgamma(n / 2.0 + 1)

        # If B is a BKZ reduced basis, this returns an estimate for the shortest vector in the lattice
        # [ B | * ]
        # [ 0 |tau]
        # under the GSA assumption, where total_vol is the volume of B, and delta is the root Hermite factor.
        # if the tau is None, the instance is homogeneous, and we omit the final row/column.
        def svp_gaussian_heuristic_gsa(i, tau):
            if tau is None:
                n = d - i
                log_vol = 2 * log_projected_vol(i)
            else:
                n = d - i + 1
                log_vol = 2 * log_projected_vol(i) + 2 * log(tau)
            log_gh = 1.0 / n * (log_vol - 2 * ball_log_vol(n))
            return log_gh

        if d > 4096:
            # chosen since RC.ADPS16(1754, 1754).log(2.) = 512.168000000000
            min_i = d - 1754
        else:
            min_i = 0

        if is_homogeneous:
            tau = None
            for i in range(min_i, d):
                if svp_gaussian_heuristic_gsa(i, tau) < log(D.stddev**2 * (d - i)):
                    return ZZ(d - (i - 1))
            return ZZ(2)
        else:
            # we look for the largest i such that (pi_i(e), tau) is shortest in the embedding lattice
            # [pi_i(B) | * ]
            # [   0    |tau]
            tau = D.stddev
            for i in range(min_i, d):
                if svp_gaussian_heuristic_gsa(i, tau) < log(D.stddev**2 * (d - i) + tau ** 2):
                    return ZZ(d - (i - 1) + 1)
            return ZZ(2)

    @classmethod
    @cached_function
    def beta_params(
        cls,
        beta: int,
        params: LWEParameters,
        zeta: int = 0,
        babai=False,
        mitm=False,
        m: int = oo,
        d: int = None,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=5,
    ):
        '''
        The costs in a Primal Hybrid attack when we run BKZ-β on the lattice basis.

        :param beta: blocksize.
        :param params: LWE parameters.
        :param zeta: guessing dimension.
        :param babai: Insist on Babai's algorithm for finding close vectors.
        :param mitm: Simulate MITM approach (√ of search space).
        :param m: number of LWE samples.
        :param d: rank of lattice to BKZ reduce. If None, we calculate the optimal dimension.

        We return a dictionary of the following values:

        - bkz_cost: the cost of BKZ-β according to the cost model
        - d: the lattice rank.
        - svp_cost: the cost of the CVP call when we use this β, according to babai=True/False.
        - eta: the projection dimension.
        - babai_probability: the probability the Babai lift in the CVP subroutine succeeds.
        - mitm_probability: the probability a mitm speedup succeeds. If mitm=False, returns 1.
        '''
        if m - zeta < beta:
            # cannot BKZ-β on a basis of dimension < β
            return {"bkz_cost": Cost(rop=oo)}

        if d is not None and d < beta:
            # cannot BKZ-β on a basis of dimension < β
            return {"bkz_cost": Cost(rop=oo)}

        simulator = simulator_normalize(red_shape_model)
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)

        if d is None:
            delta = deltaf(beta)
            d = max(beta, min(ceil(sqrt((params.n - zeta) * log(params.q / xi) / log(delta))), m - zeta))

        # 1. Simulate BKZ-β
        # We simulate BKZ-β on the dxd basis B_BKZ:
        # [q I_m |  A_{n - zeta}  ]
        # [  0   | xi I_{n - zeta}]
        # r holds the simulated squared GSO norms after BKZ-β
        r = simulator(d, params.n - zeta, params.q, beta, xi=xi, tau=False, dual=True)
        bkz_cost = costf(red_cost_model, beta, d)

        # 2. Required SVP dimension η + 1
        # We select η such that (pi_{d - η + 1}(e | s_{n - zeta}), tau) is the shortest vector in
        # [pi(B_BKZ) | t ]
        # [    0     |tau]
        if babai:
            eta = 2
            svp_cost = PrimalHybrid.babai_cost(d)
        else:
            # we scaled the lattice so that χ_e is what we want
            if red_shape_model == "gsa":
                log_vol = RR((d - (params.n - zeta)) * log(params.q) + (params.n - zeta) * log(xi))
                log_delta = RR(log(deltaf(beta)))
                svp_dim = PrimalHybrid.svp_dimension_gsa(d, log_vol, log_delta, params.Xe, params._homogeneous)
            else:
                svp_dim = PrimalHybrid.svp_dimension(r, params.Xe, is_homogeneous=params._homogeneous)
            eta = svp_dim if params._homogeneous else svp_dim - 1
            if eta > d:
                # Lattice reduction was not strong enough to "reveal" the LWE solution.
                # A larger `beta` should perhaps be attempted.
                return {"svp_cost": Cost(rop=oo)}
            # we make one svp call on a lattice of rank eta + 1
            svp_cost = costf(red_cost_model, svp_dim, svp_dim)
            # when η ≪ β, lifting may be a bigger cost
            svp_cost["rop"] += PrimalHybrid.babai_cost(d - eta)["rop"]

        if babai:
            babai_probability = prob_babai(r, sqrt(d) * params.Xe.stddev)
        else:
            babai_probability = prob_babai(r[:d-eta], sqrt(d - eta) * params.Xe.stddev)

        if mitm and zeta > 0:
            if babai:
                mitm_probability = mitm_babai_probability(r, params.Xe.stddev)
            else:
                # TODO: the probability in this case needs to be analysed
                mitm_probability = 1
        else:
            mitm_probability = 1

        return {"bkz_cost": bkz_cost,
                "d": d,
                "svp_cost": svp_cost,
                "eta": eta,
                "babai_probability": babai_probability,
                "mitm_probability": mitm_probability}

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
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        search_space=None,
        hit_probability=None,
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
        :param search_space: the size of the search space in a primal hybrid attack.
        :param hit_probability: the probability this search space hits the secret.

        .. note :: This is the lowest level function that runs no optimization, it merely reports
           costs.

        """
        beta_params = PrimalHybrid.beta_params(beta=beta,
                                               params=params,
                                               zeta=zeta,
                                               babai=babai,
                                               mitm=mitm, m=m,
                                               d=d,
                                               red_shape_model=red_shape_model,
                                               red_cost_model=red_cost_model)
        if len(beta_params.keys()) == 1:
            # this beta is not sufficient to reveal the error for these params,
            # either due to insufficient samples or projection dim > d.
            return Cost(rop=oo)

        # Search
        # MITM or no MITM
        # TODO: this is rather clumsy as a model
        def ssf(x):
            if mitm:
                return RR(sqrt(x))
            else:
                return x

        # if no search_space provided, we determine the optimal one recursively
        if search_space is None or hit_probability is None:
            f = partial(
                PrimalHybrid.cost,
                beta=beta,
                params=params,
                zeta=zeta,
                babai=babai,
                mitm=mitm,
                m=m,
                d=beta_params["d"],
                red_shape_model=red_shape_model,
                red_cost_model=red_cost_model,
            )
            min_hw = max(0, zeta - params.n + params.Xs.hamming_weight)
            max_hw = min(params.Xs.hamming_weight, zeta)
            cost = Cost(rop=oo)
            for hw in range(min_hw, max_hw + 1):
                search_space, hit_probability = guessing_set_and_hit_probability(zeta, params.Xs, hw)
                new_cost = f(search_space=search_space, hit_probability=hit_probability)
                if new_cost["rop"] > cost["rop"]:
                    # cost has started increasing, time to stop
                    return cost
                cost = new_cost
            return cost

        else:
            # we have the search_space and hit probability
            svp_cost = beta_params["svp_cost"].repeat(ssf(search_space))
            probability = hit_probability

        probability *= beta_params["babai_probability"]
        probability *= beta_params["mitm_probability"]

        bkz_cost = beta_params["bkz_cost"]
        ret = Cost()
        ret["rop"] = bkz_cost["rop"] + svp_cost["rop"]
        ret["red"] = bkz_cost["rop"]
        ret["svp"] = svp_cost["rop"]
        ret["beta"] = beta
        ret["eta"] = beta_params["eta"]
        ret["zeta"] = zeta
        ret["|S|"] = search_space
        ret["d"] = beta_params["d"]
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

        # Repeat whole experiment ~1/prob times
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
        red_shape_model=red_shape_model_default,
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
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            m=m,
            **kwds,
        )

        if baseline_cost["rop"] == oo:
            # these parameters mean usvp does not succeed for any beta < max_beta_global,
            # so we search over the full beta range
            max_beta = max_beta_global
        else:
            max_beta = baseline_cost["beta"]

        # step 1. optimize β. If zeta > 0, optimize search space.
        # the cost curve with beta is non-smooth, with sudden jumps due to the search space changing. We eliminate
        # these jumps by instead fixing the search space, and finding the best attack for that search space.
        # we then loop over increasing search spaces to find the best overall attack.
        # our search space is formed of all zeta length strings of weight at most hw for some hw.
        # the smallest admissable hw
        min_hw = max(0, zeta - params.Xs.n + params.Xs.hamming_weight)
        # the largest admissable hw
        max_hw = min(zeta, params.Xs.hamming_weight)
        cost = Cost(rop=oo)
        for hw in range(min_hw, max_hw + 1):
            search_space, hit_probability = guessing_set_and_hit_probability(zeta, params.Xs, hw)
            precision = 2
            with local_minimum(40, max_beta + precision, precision=precision, log_level=log_level + 1) as it:
                for beta in it:
                    it.update(f(beta, search_space=search_space, hit_probability=hit_probability))
                for beta in it.neighborhood:
                    it.update(f(beta, search_space=search_space, hit_probability=hit_probability))
                new_cost = it.y
            if new_cost["rop"] > cost["rop"]:
                # cost has started increasing, time to stop
                break
            else:
                cost = new_cost
        Logging.log("bdd", log_level, f"H1: {cost!r}")

        # step 2. optimize d
        if cost and cost.get("tag", "XXX") != "usvp" and optimize_d:
            with local_minimum(
                params.n - zeta, cost["d"] + 1, log_level=log_level + 1
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
            rop: ≈2^87.9, red: ≈2^87.2, svp: ≈2^86.5, β: 127, η: 18, ζ: 303, |S|: ≈2^45.9, d: 409, prob: ≈2^-19.4, ↻...

            >>> LWE.primal_hybrid(params, mitm=False, babai=True)
            rop: ≈2^88.2, red: ≈2^87.3, svp: ≈2^87.1, β: 119, η: 2, ζ: 304, |S|: ≈2^45.9, d: 401, prob: ≈2^-21.6, ↻...

            >>> LWE.primal_hybrid(params, mitm=True, babai=False)
            rop: ≈2^72.1, red: ≈2^70.8, svp: ≈2^71.4, β: 106, η: 18, ζ: 321, |S|: ≈2^82.8, d: 373, prob: 0.002, ↻...

            >>> LWE.primal_hybrid(params, mitm=True, babai=True)
            rop: ≈2^85.4, red: ≈2^83.6, svp: ≈2^84.9, β: 111, η: 2, ζ: 366, |S|: ≈2^90.9, d: 330, prob: ≈2^-20.5, ↻...

        TESTS:

        We test a trivial instance::

            >>> params = LWE.Parameters(2**10, 2**100, ND.DiscreteGaussian(3.19), ND.DiscreteGaussian(3.19))
            >>> LWE.primal_bdd(params)
            rop: ≈2^43.6, red: ≈2^43.6, svp: ≈2^21.9, β: 40, η: 2, d: 1511, tag: bdd

        We also test a LWE instance with a large error (coming from issue #106)::

            >>> LWE.primal_bdd(LWE.Parameters(n=256, q=12289, Xs=ND.UniformMod(2), Xe=ND.UniformMod(1024)))
            rop: ≈2^115.4, red: ≈2^76.7, svp: ≈2^115.4, β: 170, η: 336, d: 336, tag: bdd

            >>> LWE.primal_bdd(LWE.Parameters(n=700, q=2**64, Xs=ND.UniformMod(2), Xe=ND.UniformMod(2**59)))
            rop: ≈2^249.8, red: ≈2^225.4, svp: ≈2^249.8, β: 702, η: 818, d: 1023, tag: bdd

            # small n example (Issue #182)
            >>> LWE.primal_bdd(LWE.Parameters(n=8, q = 2**128, Xs=ND.UniformMod(2**128), Xe=ND.UniformMod(2**(126))))
            rop: ≈2^185.5, red: ≈2^184.4, svp: ≈2^184.7, β: 585, η: 585, d: 585, tag: bdd
        """

        if zeta == 0:
            tag = "bdd"
        else:
            tag = "hybrid"

        params = LWEParameters.normalize(params)

        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

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
            # primal_hybrid cost is generally parabolic with zeta.
            # We find a range [min_zeta, max_zeta] such that cost is finite over the entire interval.

            # we search for min_zeta such that cost(min_zeta) is finite, but cost(min_zeta - 1) is infinite.
            cost_min_zeta = f(zeta=0, optimize_d=False, **kwds)
            if cost_min_zeta["rop"] < oo:
                min_zeta = 0
            else:
                min_zeta_lower = 0
                min_zeta_upper = params.n - 1
                min_zeta = (min_zeta_upper + min_zeta_lower) // 2
                while min_zeta_upper - min_zeta_lower > 1:
                    cost_min_zeta = f(min_zeta, optimize_d=False, **kwds)
                    if cost_min_zeta["rop"] < oo:
                        min_zeta_upper = min_zeta
                    else:
                        min_zeta_lower = min_zeta
                    min_zeta = (min_zeta_upper + min_zeta_lower) // 2

            # we search for max_zeta such that cost(max_zeta) is finite, but cost(max_zeta + 1) is infinite.
            cost_max_zeta = f(zeta=params.n-1, optimize_d=False, **kwds)
            if cost_max_zeta["rop"] < oo:
                max_zeta = params.n - 1
            else:
                max_zeta_lower = 0
                max_zeta_upper = params.n - 1
                max_zeta = (max_zeta_upper + max_zeta_lower) // 2
                while max_zeta_upper - max_zeta_lower > 1:
                    cost_max_zeta = f(max_zeta, optimize_d=False, **kwds)
                    if cost_max_zeta["rop"] < oo:
                        max_zeta_lower = max_zeta
                    else:
                        max_zeta_upper = max_zeta
                    max_zeta = (max_zeta_upper + max_zeta_lower) // 2

            ret = minimize_scalar(lambda x: log(f(zeta=round(x), optimize_d=False,
                                                  **kwds)["rop"]), bounds=(min_zeta, max_zeta), method="bounded")

            zeta = int(ret.x)
            cost = f(zeta=zeta, optimize_d=False, **kwds)
            # check a small neighborhood of this zeta
            precision = 3
            for zeta in range(max(0, zeta - precision), min(params.n, zeta + precision) + 1):
                cost = min(cost, f(zeta=zeta, optimize_d=False, **kwds))
            # minimize_scalar fits to a parabola. This can cause this search to miss minima at extrema
            cost = min(cost, cost_min_zeta, cost_max_zeta)

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
