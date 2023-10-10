# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from functools import partial

from sage.all import oo, ceil, sqrt, log, RR, cached_function, exp, pi, floor
from math import lgamma
from scipy.special import digamma
import numpy as np
from .reduction import delta as deltaf
from .reduction import cost as costf
from .util import local_minimum, zeta_precomputed, zeta_prime_precomputed, gh_constant
from .cost import Cost
from .lwe_primal import PrimalUSVP, PrimalHybrid
from .ntru_parameters import NTRUParameters
from .simulator import normalize as simulator_normalize
from .prob import conditional_chi_squared, chisquared_table
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .conf import red_simulator as red_simulator_default
from .conf import max_n_cache


class PrimalDSD():
    """
    Estimate cost of solving (overstretched) NTRU via dense sublattice discovery
    """
    @staticmethod
    @cached_function
    def ball_log_vol(n):
        return RR((n/2.)) * RR(log(pi)) - RR(lgamma(n/2. + 1))

    @staticmethod
    def log_gh(d, logvol=0):
        if d < 49:
            return RR(gh_constant[d]) + RR(logvol)/d

        return RR(1./d) * RR(logvol - PrimalDSD.ball_log_vol(d))

    @staticmethod
    def DSL_logvol_matrix(n, sigmasq):
        total = n*(RR(log(sigmasq))+RR(log(2.))+RR(digamma(n)))/2.
        proj_loss = np.sum([(digamma((2*n-i)/2.)-digamma(n)) for i in range(n)])/2.
        return total+proj_loss

    @staticmethod
    def DSL_logvol_circulant(n, sigmasq):
        lambda0 = RR((np.log(2)-np.euler_gamma+np.log(n)+np.log(sigmasq))/2.)
        lambdai = (n-1)*(1-np.euler_gamma+np.log(n)+np.log(sigmasq))/2.
        return lambda0+lambdai

    @staticmethod
    def DSL_logvol_circulant_fixed(n, R):
        lambda0 = (-np.euler_gamma+np.log(R))/2.
        lambdai = (n-1)*(1-np.euler_gamma+np.log(R)-np.log(2))/2.
        return lambda0+lambdai

    @staticmethod
    @cached_function
    def DSL_logvol(n, sigmasq, ntru="circulant"):
        if ntru=="matrix":
            return PrimalDSD.DSL_logvol_matrix(n, sigmasq)
        if ntru=="circulant":
            return PrimalDSD.DSL_logvol_circulant(n, sigmasq)
        if ntru=="fixed":
            return PrimalDSD.DSL_logvol_circulant_fixed(n, sigmasq)
        print("non implemented ntru type")

    @staticmethod
    @cached_function
    def proj_logloss(d, k):
        # log loss of length when projecting out k dimension out of d
        return (RR(digamma((d-k)/2.))-RR(digamma(d/2.)))/2.

    @staticmethod
    def DSLI_vols(dsl_logvol, FL_shape):
        n = len(FL_shape)//2
        vols = (2*n+1)*[None]

        dsl_dim = n
        vols[2*n] = dsl_logvol

        # Going to a intersection of dimension s
        for s in range(2*n-1, n, -1):
            # Negate cause it's a dual vector really
            x = - FL_shape[s]
            x += PrimalDSD.proj_logloss(s+1, n)
            x += zeta_prime_precomputed[dsl_dim]/zeta_precomputed[dsl_dim]  # primitivity
            dsl_logvol += x
            vols[s] = dsl_logvol
            dsl_dim -= 1

        assert dsl_dim == 1
        assert s == n+1

        return vols

    @staticmethod
    @cached_function
    def prob_dsd(
        beta: int,
        params: NTRUParameters,
        simulator,
        m: int = oo,
        tau=None,
        d=None,
        dsl_logvol=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
    ):

        if d is None:
            d = m

        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)
        if dsl_logvol is None:
            dsl_logvol = PrimalDSD.DSL_logvol(params.n, params.Xs.stddev**2, ntru=params.ntru_type)

        B_shape = [log(r_)/2 for r_ in simulator(d, params.n, params.q, beta, xi=xi, tau=tau)]

        dsli_vols = PrimalDSD.DSLI_vols(dsl_logvol, B_shape)

        prob_all_not = RR(1.)
        prob_pos = np.zeros(2*params.n, dtype='double')
        for i in range(1, params.n+1):
            s = params.n + i

            dslv_len = PrimalDSD.log_gh(i, dsli_vols[s])
            sigma_sq = exp(2*dslv_len)/s

            if sigma_sq > 10**10:
                prob_pos[s-beta] = 0.
                continue

            norm_threshold = exp(2*(B_shape[s-beta]))/sigma_sq
            proba_one = chisquared_table[beta].cum_distribution_function(norm_threshold)

            if proba_one <= 10e-8:
                continue

            # account for pulling back probability if beta small
            if beta <= 20:
                for j in range(2, int(s/beta+1)):
                    if proba_one < 10**(-6):
                        proba_one = 0.
                        break
                    ind = s - j*(beta-1)-1
                    norm_bt = exp(2*B_shape[ind])/sigma_sq
                    norm_b2 = exp(2*B_shape[ind+beta-1])/sigma_sq
                    proba_one *= conditional_chi_squared(beta-1, s-ind-(beta-1), norm_bt, norm_b2)

            prob_pos[s-beta] = proba_one
            prob_all_not *= max(1.-proba_one, 0.)
            Logging.log("dsd", log_level+1, f"Pr[dsd, {beta}] = {prob_all_not}")

        return RR(1.-prob_all_not), prob_pos

    def __call__(
        self,
        params: NTRUParameters,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        log_level=1,
        **kwds,
    ):
        """
        Estimate cost of solving (overstretched) NTRU using the Dense sublattice.
        Code is adapted from Léo Ducas and Wessel van Woerden.
        See https://github.com/WvanWoerden/NTRUFatigue for original source

        :param params: NTRU parameters.
        :param red_cost_model: How to cost lattice reduction.
        :param red_shape_model: How to model the shape of a reduced basis.
        :return: A cost dictionary.

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction.
        - ``δ``: Root-Hermite factor targeted by lattice reduction.
        - ``β``: BKZ block size.
        - ``d``: Lattice dimension.

        EXAMPLE::

            >>> from estimator import *
            >>> NTRU.primal_dsd(schemes.NTRUHPS2048509Enc)
            rop: ≈2^157.1, red: ≈2^157.1, δ: 1.003645, β: 453, d: 1016, tag: dsd

            NOTE: Non-overstretched parameters (where the probability of Dense sublattice
            discovery is 0) will return β = d.

            >>> params = NTRU.Parameters(n=113, q=512, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            >>> NTRU.primal_dsd(params, red_shape_model="zgsa")
            rop: ≈2^41.3, red: ≈2^41.3, δ: 1.012468, β: 42, d: 226, tag: dsd

            >>> NTRU.primal_dsd(params, red_shape_model="cn11")
            rop: ≈2^41.2, red: ≈2^41.2, δ: 1.012468, β: 41, d: 226, tag: dsd

            >>> NTRU.primal_dsd(params, red_shape_model=Simulator.CN11)
            rop: ≈2^41.2, red: ≈2^41.2, δ: 1.012468, β: 41, d: 226, tag: dsd

        The success condition was formulated in [EC:KF17] and further refined in [AC:DvW21]
        """
        params = NTRUParameters.normalize(params)
        # allow for a larger embedding lattice dimension: Bai and Galbraith
        m = params.m + params.n if params.Xs <= params.Xe else params.m

        try:
            red_shape_model = simulator_normalize(red_shape_model)
        except ValueError:
            pass

        if params.n > max_n_cache:
            raise ValueError("Please increase the hardcoded value of max_n_cache to run the predictor for such large n")

        remaining_proba = RR(1.)
        average_beta = RR(0.)
        total_DSD_prob = RR(0.)
        DSD_prob = RR(0.)
        prob_pos_total = np.zeros(2*params.n, dtype='double')

        for beta in range(2, params.n):
            tours = floor(params.n**2 / beta**2)+3

            DSD_prob, DSD_prob_pos = self.prob_dsd(beta, params, red_shape_model, m=m,
                                                   red_cost_model=red_cost_model, log_level=log_level)

            if DSD_prob > 10e-8:
                for t in range(tours):
                    for i in range(2*params.n):
                        prob_pos = DSD_prob_pos[i]
                        average_beta += RR(beta) * remaining_proba * prob_pos
                        prob_pos_total[i] += remaining_proba * prob_pos
                        total_DSD_prob += remaining_proba * prob_pos
                        remaining_proba *= (1.-prob_pos)

                Logging.log("dsd", log_level+1, "β= %d,\t pr=%.4e, \t rem-pr=%.4e"%(beta, DSD_prob, remaining_proba))
            if remaining_proba < 0.001:
                average_beta += beta * remaining_proba
                break

        if not average_beta:
            average_beta = m
            predicate = False

        else:
            predicate = True

        cost = costf(red_cost_model, RR(average_beta), m, predicate=predicate)
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
