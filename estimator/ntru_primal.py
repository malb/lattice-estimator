# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from sage.all import oo, log, RR, cached_function, exp, pi, floor, euler_gamma
from math import lgamma
from scipy.special import digamma
from .reduction import cost as costf
from .util import zeta_precomputed, zeta_prime_precomputed, gh_constant
from .lwe_primal import PrimalUSVP, PrimalHybrid
from .ntru_parameters import NTRUParameters
from .simulator import normalize as simulator_normalize
from .prob import conditional_chi_squared, chisquared_table
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
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
        proj_loss = sum([(digamma((2*n-i)/2.)-digamma(n)) for i in range(n)])/2.
        return total+proj_loss

    @staticmethod
    def DSL_logvol_circulant(n, sigmasq):
        lambda0 = RR((log(2) - euler_gamma + log(n) + log(sigmasq))/2.)
        lambdai = RR((n - 1)*(1 - euler_gamma + log(n) + log(sigmasq))/2.)
        return lambda0+lambdai

    @staticmethod
    def DSL_logvol_circulant_fixed(n, R):
        lambda0 = RR((-euler_gamma + log(R))/2.)
        lambdai = RR((n - 1)*(1 - euler_gamma + log(R) - log(2))/2.)
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

        raise ValueError(f"NTRU type: {ntru} is not supported.")

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
        prob_pos = (2*params.n)*[RR(0)]
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
            >>> NTRU.primal_dsd(schemes.NTRUHRSS701Enc)
            rop: ≈2^190.2, red: ≈2^190.2, δ: 1.003095, β: 571, d: 1400, tag: dsd

            >>> params = NTRU.Parameters(n=113, q=512, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            >>> NTRU.primal_dsd(params, red_shape_model="zgsa")
            rop: ≈2^41.3, red: ≈2^41.3, δ: 1.012468, β: 42, d: 226, tag: dsd

            >>> NTRU.primal_dsd(params, red_shape_model="cn11")
            rop: ≈2^41.2, red: ≈2^41.2, δ: 1.012468, β: 41, d: 226, tag: dsd

            >>> NTRU.primal_dsd(params, red_shape_model=Simulator.CN11)
            rop: ≈2^41.2, red: ≈2^41.2, δ: 1.012468, β: 41, d: 226, tag: dsd

        The success condition was formulated in [EC:KirFou17]_ and further refined in [AC:DucWoe21]_

        .. note :: Non-overstretched parameters (where the probability of Dense sublattice
           discovery is 0) will return β = d.
        """
        if params.Xs.stddev != params.Xe.stddev:
            raise NotImplementedError("Dense sublattice attack not supported for Xs != Xe")

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
        prob_pos_total = (2*params.n)*[RR(0.)]

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

    def __call__(
        self,
        params: NTRUParameters,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        optimize_d=True,
        log_level=1,
        **kwds,
    ):
        """
        Estimate cost of solving NTRU via uSVP reduction.

        :param params: NTRU parameters.
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
            >>> NTRU.primal_usvp(schemes.NTRUHPS2048509Enc)
            rop: ≈2^134.6, red: ≈2^134.6, δ: 1.004179, β: 373, d: 929, tag: usvp

            >>> params = NTRU.Parameters(n=200, q=127, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
            >>> NTRU.primal_usvp(params, red_shape_model="cn11")
            rop: ≈2^87.2, red: ≈2^87.2, δ: 1.006132, β: 208, d: 374, tag: usvp

            >>> NTRU.primal_usvp(params, red_shape_model=Simulator.CN11)
            rop: ≈2^87.2, red: ≈2^87.2, δ: 1.006132, β: 208, d: 374, tag: usvp

            >>> NTRU.primal_usvp(params, red_shape_model=Simulator.CN11, optimize_d=False)
            rop: ≈2^87.4, red: ≈2^87.4, δ: 1.006132, β: 208, d: 399, tag: usvp

        The success condition was formulated in [USENIX:ADPS16]_ and studied/verified in
        [AC:AGVW17]_, [C:DDGR20]_, [PKC:PosVir21]_. The treatment of small secrets is from
        [ACISP:BaiGal14]_.

        """
        return super().__call__(params,
                                red_cost_model=red_cost_model,
                                red_shape_model=red_shape_model,
                                optimize_d=optimize_d,
                                log_level=log_level,
                                **kwds)


primal_usvp = NTRUPrimalUSVP()


class NTRUPrimalHybrid(PrimalHybrid):

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
        """
        Estimate the cost of the hybrid attack and its variants.

        :param params: NTRU parameters.
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
            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = False, babai = False)
            rop: ≈2^87.8, red: ≈2^87.0, svp: ≈2^86.6, β: 116, η: 21, ζ: 302, |S|: ≈2^39.2, d: 372, prob: ≈2^-22.3, ↻...

            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = False, babai = True)
            rop: ≈2^88.0, red: ≈2^87.4, svp: ≈2^86.4, β: 98, η: 2, ζ: 318, |S|: ≈2^39.6, d: 328, prob: ≈2^-27.9, ↻: ...

            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = True, babai = False)
            rop: ≈2^80.1, red: ≈2^79.7, svp: ≈2^78.3, β: 170, η: 22, ζ: 254, |S|: ≈2^103.7, d: 495, prob: 0.708, ↻: ...

            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = True, babai = True)
            rop: ≈2^85.1, red: ≈2^84.1, svp: ≈2^84.0, β: 105, η: 2, ζ: 363, |S|: ≈2^85.0, d: 294, prob: ≈2^-22.9, ↻:...

        TESTS:

        We test a trivial instance::

            >>> params = NTRU.Parameters(2**10, 2**100, ND.DiscreteGaussian(3.19), ND.DiscreteGaussian(3.19))
            >>> NTRU.primal_bdd(params)
            rop: ≈2^43.6, red: ≈2^43.6, svp: ≈2^35.0, β: 40, η: 46, d: 1461, tag: bdd

        """
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
