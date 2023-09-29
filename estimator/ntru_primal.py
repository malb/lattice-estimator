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
from .lwe_primal import PrimalUSVP, PrimalHybrid
from .ntru_parameters import NTRUParameters
from .overstretched_ntru import find_fatigue, combined_attack_prob
from .errors import NotOverstretchedError
from .simulator import normalize as simulator_normalize
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai
from .prob import mitm_babai_probability
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .conf import red_simulator as red_simulator_default

class PrimalDSD():
    """
    Estimate cost of solving (overstretched) NTRU via uSVP reduction
    """

    @staticmethod
    @cached_function
    def cost_simulator(
        beta: int,
        params: NTRUParameters,
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

        # TODO: this predicate isn't quite right for DSD estimate.
        # TODO: In fact, we could replace cost simulator directly by the DSD_attack_prob function
        # TODO: Remove the overstretched NTRU file.

        # r = simulator(d=d, n=params.n, q=params.q, beta=beta, xi=xi, tau=tau)
        # lhs = params.Xe.stddev**2 * (beta - 1)
        # predicate = r[d - beta] > lhs
        if not beta:
            beta = d
            predicate = False
        
        else: 
            predicate = True

        return costf(red_cost_model, beta, d, predicate=predicate)

    def __call__(
        self,
        params: NTRUParameters,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        log_level=1,
        **kwds,      
    ):
        """
        TODO Add parameter description
        """
        params = NTRUParameters.normalize(params)

        try:
            red_shape_model = simulator_normalize(red_shape_model)
        except ValueError:
            pass
        
        # Define Basis Shape Function for NTRU fatigue estimator
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        zshapef = partial(red_shape_model,
                    d=(params.n + params.m),
                    xi=xi, tau=None)
        
        # Use NTRU fatigue Estimator for DSD event hardness
        beta, sk_prob, dsd_prob, _ = combined_attack_prob(params.q, params.n, sk_variance=params.Xs.stddev**2, ntru=params.ntru_type, verbose=False,
                                             only="DSD", zshapef=zshapef)

        # print(beta, sk_prob, dsd_prob)
        cost = self.cost_simulator(int(beta), params, red_shape_model, m=params.m, 
                                   d=(params.m + params.n), red_cost_model=red_cost_model)
        cost["tag"] = "dsd"
        cost["problem"] = params
        return cost.sanity_check()

    __name__ = "primal_dsd"

primal_dsd = PrimalDSD()

class NTRUPrimalUSVP(PrimalUSVP):
    
    @staticmethod
    def _solve_for_d(params, m, beta, tau, xi):
        """
        Find smallest d ∈ [n,m] to satisfy uSVP condition.

        If no such d exists, return the upper bound m.
        :tau: ignored due to NTRU being homogeneous
        """
        # Find the smallest d ∈ [n,m] s.t. a*d^2 + b*d + c >= 0

        delta = deltaf(beta)
        a = -log(delta)
        C = log(params.Xe.stddev**2 * (beta - 1)) / 2.0
        b = log(delta) * (2 * beta - 1) + log(params.q) - C
        c = params.n * log(xi) - (params.n + 1) * log(params.q)
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
        params: NTRUParameters,
        m: int = oo,
        tau=None,
        d=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
    ):
        """
        GSA Cost function modified to operate on a homogeneous (i.e. b=0) instance.
        """
        delta = deltaf(beta)
        xi = NTRUPrimalUSVP._xi_factor(params.Xs, params.Xe)
        m = min(2 * ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        d = NTRUPrimalUSVP._solve_for_d(params, m, beta, tau, xi) if d is None else d
        # if d == β we assume one SVP call, otherwise poly calls. This makes the cost curve jump, so
        # we avoid it here.
        if d == beta and d < m:
            d += 1
        assert d <= m + 1

        lhs = log(sqrt(params.Xe.stddev**2 * (beta - 1)))
        rhs = RR(
            log(delta) * (2 * beta - d - 1)
            + (log(xi) * params.n + log(params.q) * (d - params.n - 1)) / d
        )

        return costf(red_cost_model, beta, d, predicate=lhs <= rhs)

    @staticmethod
    @cached_function
    def cost_simulator(
        beta: int,
        params: NTRUParameters,
        simulator,
        m: int = oo,
        tau=None,
        d=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
    ):
        """
        cost simulator function modified to operate on a homogeneous (i.e. b=0) instance.
        """
        delta = deltaf(beta)
        if d is None:
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        xi = NTRUPrimalUSVP._xi_factor(params.Xs, params.Xe)

        r = simulator(d=d, n=params.n, q=params.q, beta=beta, xi=xi, tau=tau)
        lhs = params.Xe.stddev**2 * (beta - 1)
        predicate = r[d - beta] > lhs

        return costf(red_cost_model, beta, d, predicate=predicate)
    
    def __call__(
        self,
        params: NTRUParameters,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        optimize_d=True,
        log_level=1,
        **kwds,
    ):
        # breakpoint()
        return super().__call__(params,
                         red_cost_model=red_cost_model,
                         red_shape_model=red_shape_model,
                         optimize_d=optimize_d,
                         log_level=log_level,
                         **kwds)

primal_usvp = NTRUPrimalUSVP()

class NTRUPrimalHybrid(PrimalHybrid):

    @classmethod
    def cost_zeta(
        cls,
        zeta: int,
        params: NTRUParameters,
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
        params: NTRUParameters,
        babai: bool = True,
        zeta: int = None,
        mitm: bool = True,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=1,
        **kwds,
    ):
        return super().__call__(params,
                         babai=babai,
                         zeta=zeta,
                         mitm=mitm,
                         red_shape_model=red_shape_model,
                         red_cost_model=red_cost_model,
                         log_level=log_level,
                         **kwds)
    
primal_hybrid = NTRUPrimalHybrid()

def primal_bdd(
    params: NTRUParameters,
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