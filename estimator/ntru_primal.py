# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using primal attacks.

See :ref:`LWE Primal Attacks` for an introduction what is available.

"""
from sage.all import oo, ceil, sqrt, log, RR, cached_function, exp, pi, floor, matrix, euler_gamma, binomial
from math import lgamma
from scipy.special import digamma
from .reduction import delta as deltaf
from .reduction import cost as costf
from .util import zeta_precomputed, zeta_prime_precomputed, gh_constant
from .cost import Cost
from .lwe_primal import PrimalUSVP, PrimalHybrid
from .ntru_parameters import NTRUParameters
from .simulator import normalize as simulator_normalize
from .prob import conditional_chi_squared, chisquared_table
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .prob import babai as prob_babai
from .prob import mitm_babai_probability
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
            >>> NTRU.primal_dsd(schemes.NTRUHPS2048509Enc)
            rop: ≈2^157.1, red: ≈2^157.1, δ: 1.003645, β: 453, d: 1016, tag: dsd

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

        .. note :: This function is overloaded from LWE.PrimalUSVP to account for the fact that
           NTRU is homogeneous
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

        .. note :: This function is overloaded from LWE.PrimalUSVP to account for the fact that
           NTRU is homogeneous.
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
            rop: ≈2^134.5, red: ≈2^134.5, δ: 1.004179, β: 373, d: 923, tag: usvp

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

    @staticmethod
    @cached_function
    def cost(
        beta: int,
        params: NTRUParameters,
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

        .. note :: This function is overloaded from LWE.PrimalHybrid to account for the fact that
           NTRU is homogeneous

        """
        if d is None:
            delta = deltaf(beta)
            d = min(ceil(sqrt(params.n * log(params.q) / log(delta))), m)
        d -= zeta
        xi = PrimalUSVP._xi_factor(params.Xs, params.Xe)

        # 1. Simulate BKZ-β
        # tau is None due to homogenous nature of NTRU

        r = simulator(d, params.n - zeta, params.q, beta, xi=xi, tau=None, dual=True)

        bkz_cost = costf(red_cost_model, beta, d)

        # 2. Required SVP dimension η
        if babai:
            eta = 2
            svp_cost = PrimalHybrid.babai_cost(d)
        else:
            # we scaled the lattice so that χ_e is what we want
            eta = PrimalHybrid.svp_dimension(r, params.Xe)
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
            h = ceil(len(params.Xs) * params.Xs.density)
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
                probability *= mitm_babai_probability(r, params.Xe.stddev, params.q)
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
            rop: ≈2^87.8, red: ≈2^87.0, svp: ≈2^86.6, β: 116, η: 21, ζ: 302, |S|: ≈2^39.2, d: 372, prob: ≈2^-22.3, ...

            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = False, babai = True)
            rop: ≈2^88.0, red: ≈2^87.4, svp: ≈2^86.4, β: 98, η: 2, ζ: 318, |S|: ≈2^39.6, d: 328, prob: ≈2^-27.9, ...

            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = True, babai = False)
            rop: ≈2^80.1, red: ≈2^79.7, svp: ≈2^78.3, β: 170, η: 22, ζ: 254, |S|: ≈2^103.7, d: 495, prob: 0.708, ...

            >>> NTRU.primal_hybrid(schemes.NTRUHPS2048509Enc.updated(Xs=ND.SparseTernary(508,16)),
            ... mitm = True, babai = True)
            rop: ≈2^85.1, red: ≈2^84.1, svp: ≈2^84.0, β: 105, η: 2, ζ: 363, |S|: ≈2^85.0, d: 294, prob: ≈2^-22.9, ...

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
