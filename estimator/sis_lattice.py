# -*- coding: utf-8 -*-
"""
Estimate cost of solving SIS using lattice reduction attacks.

See :ref:`SIS Lattice Attacks` for an introduction what is available.

"""
from functools import partial

from sage.all import oo, sqrt, log, RR, floor, cached_function
from .reduction import beta as betaf
from .reduction import cost as costf
from .util import local_minimum
from .cost import Cost
from .sis_parameters import SISParameters
from .simulator import normalize as simulator_normalize
from .prob import gaussian_cdf
from .prob import amplify as prob_amplify
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .conf import red_simulator as red_simulator_default


class SISLattice:
    """
    Estimate cost of solving SIS via lattice reduction.
    """

    @staticmethod
    def _solve_for_delta_euclidean(params, d):
        # root_volume = params.q**(params.n/d)
        # delta = (params.length_bound / root_volume)**(1/(d - 1))
        root_volume = (params.n / d) * log(params.q, 2)
        log_delta = (1 / (d - 1)) * (log(params.length_bound, 2) - root_volume)
        return RR(2**log_delta)

    @staticmethod
    def _opt_sis_d(params):
        """
        Optimizes SIS dimension for the given parameters, assuming the optimal
        d \approx sqrt(n*log(q)/log(delta))
        """
        log_delta = log(params.length_bound, 2) ** 2 / (4 * params.n * log(params.q, 2))
        d = sqrt(params.n * log(params.q, 2) / log_delta)
        return d

    @staticmethod
    @cached_function
    def cost_euclidean(
        params: SISParameters,
        d=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
        **kwds,
    ):
        # Check for triviality
        if params.length_bound >= params.q:
            raise ValueError("SIS trivially easy. Please set norm bound < q.")

        if d is None:
            d = min(floor(SISLattice._opt_sis_d(params)), params.m)

        # First solve for root hermite factor
        delta = SISLattice._solve_for_delta_euclidean(params, d)
        # Then derive beta from the cost model(s)
        if delta >= 1 and betaf(delta) <= d:
            beta = betaf(delta)
            reduction_possible = True

        else:
            beta = d
            reduction_possible = False

        lb = min(RR(sqrt(params.n * log(params.q))), RR(sqrt(d) * params.q ** (params.n / d)))
        return costf(
            red_cost_model, beta, d, predicate=params.length_bound > lb and reduction_possible
        )

    @staticmethod
    @cached_function
    def cost_infinity(
        beta: int,
        params: SISParameters,
        simulator,
        zeta: int = 0,
        success_probability: float = 0.99,
        d=None,
        red_cost_model=red_cost_model_default,
        log_level=None,
        **kwds,
    ):
        """
        Computes the cost of the attack on SIS using an infinity norm bound.

        :param params: SIS parameters
        :param beta: Block size used to produce short vectors for reduction
        :param simulator: Basis profile simulator
        :param zeta: Number of coefficients to set to 0 (ignore)
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction

        .. note :: This function assumes that the instance is normalized. It runs no optimization,
            it merely reports costs.

        """
        if params.length_bound >= params.q:
            raise ValueError("SIS trivially easy. Please set norm bound < q.")

        if d is None:
            d = params.m

        # Calculate the basis shape to aid in both styles of analysis
        d_ = d - zeta

        if d_ < beta:
            return Cost(rop=oo, mem=oo)

        r = simulator(d=d_, n=d_ - params.n, q=params.q, beta=beta, xi=1, tau=False)

        # Cost the sampling of short vectors.
        rho, cost_red, N, sieve_dim = red_cost_model.short_vectors(beta, d_)
        bkz_cost = costf(red_cost_model, beta, d_)

        if RR(sqrt(d)) * params.length_bound <= params.q:  # Non-dilithium style analysis
            # Calculate expected vector length using approximation factor on the shortest vector from BKZ
            vector_length = rho * sqrt(r[0])
            # Find probability that all coordinates meet norm bound
            sigma = vector_length / sqrt(d_)
            log_trial_prob = RR(d_ * log(1 - 2 * gaussian_cdf(0, sigma, -params.length_bound), 2))

        else:  # Dilithium style analysis
            # Find first non-q-vector in r
            if abs(sqrt(r[0]) - params.q) < 1e-8:  # q-vectors exist
                idx_start = next(i for i, r_ in enumerate(r) if r_ < r[0])

            else:
                idx_start = 0

            if abs(r[-1] - 1) < 1e-8:  # 1-vectors exist
                # Find first 1 length graham-schmidt vector in r (Zone III)
                idx_end = next((i - 1 for i, r_ in enumerate(r) if sqrt(r_) <= 1 + 1e-8), d_ - 1)

            else:
                idx_end = d_ - 1

            vector_length = sqrt(r[idx_start])
            gaussian_coords = max(idx_end - idx_start + 1, sieve_dim)
            sigma = vector_length / sqrt(gaussian_coords)

            log_trial_prob = RR(
                log(1 - 2 * gaussian_cdf(0, sigma, -params.length_bound), 2) * (gaussian_coords)
            )
            log_trial_prob += RR(log((2 * params.length_bound + 1) / params.q, 2) * (idx_start))

        probability = 2 ** min(
            0, log_trial_prob + RR(log(N, 2))
        )  # expected number of solutions (max 1)
        ret = Cost()
        ret["rop"] = cost_red
        ret["red"] = bkz_cost["rop"]
        ret["sieve"] = max(cost_red - bkz_cost["rop"], 1e-100)  # Ensuring non-zero cost here
        ret["beta"] = beta
        ret["eta"] = sieve_dim
        ret["zeta"] = zeta
        ret["d"] = d_
        ret["prob"] = probability

        ret.register_impermanent(
            rop=True,
            red=True,
            sieve=True,
            eta=False,
            zeta=False,
            prob=False,
        )
        # 4. Repeat whole experiment ~1/prob times
        if probability and not RR(probability).is_NaN():
            ret = ret.repeat(
                prob_amplify(success_probability, probability),
            )
        else:
            return Cost(rop=oo)

        return ret

    @classmethod
    def cost_zeta(
        cls,
        zeta: int,
        params: SISParameters,
        ignore_qary: bool = False,
        red_shape_model=red_simulator_default,
        red_cost_model=red_cost_model_default,
        d=None,
        log_level=5,
        **kwds,
    ):
        """
        This function optimizes costs for a fixed number of coordinates to 'ignore', denoted ζ.
        Ignored coordinates are set to 0 in the final SIS solution, so the dimension of the
        instance is treated as d-ζ.
        """
        # step 0. establish baseline cost using worst case euclidean norm estimate
        params_baseline = params.updated(norm=2)
        baseline_cost = lattice(
            params_baseline,
            ignore_qary=ignore_qary,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            log_level=log_level + 1,
            **kwds,
        )

        Logging.log("sis_infinity", log_level, f"H0: {repr(baseline_cost)}")

        f = partial(
            cls.cost_infinity,
            params=params,
            zeta=zeta,
            ignore_qary=ignore_qary,
            simulator=red_shape_model,
            red_cost_model=red_cost_model,
            d=d,
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

        Logging.log("sis_infinity", log_level, f"H1: {cost!r}")
        if cost is None:
            return Cost(rop=oo)
        return cost

    def __call__(
        self,
        params: SISParameters,
        zeta: int = None,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=1,
        **kwds,
    ):
        """
        Estimate the cost of attacking SIS using lattice reduction

        :param params: SIS parameters.
        :param zeta: Number of coefficients to set to 0 (ignore)
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction.
        - ``δ``: Root-Hermite factor targeted by lattice reduction.
        - ``β``: BKZ block size.
        - ``η``: Dimension of the final Sieving call to generate short vectors.
        - ``ζ``: Number of ignored coordinates.
        - ``|S|``: Guessing search space.
        - ``prob``: Probability of success in guessing.
        - ``repeat``: How often to repeat the attack.
        - ``d``: Lattice dimension.

        EXAMPLES::

            >>> from estimator import *
            >>> SIS.lattice(schemes.Dilithium2_MSIS_WkUnf)
            rop: ≈2^152.2, red: ≈2^151.3, sieve: ≈2^151.1, β: 427, η: 433, ζ: 0, d: 2304, prob: 1, ↻: 1, tag: infinity

            >>> SIS.lattice(schemes.Dilithium2_MSIS_WkUnf, red_shape_model="lgsa")
            rop: ≈2^151.3, red: ≈2^150.2, sieve: ≈2^150.5, β: 423, η: 431, ζ: 0, d: 2304, prob: 1, ↻: 1, tag: infinity

            >>> params = SIS.Parameters(n=113, q=2048, length_bound=512, norm=2)
            >>> SIS.lattice(params)
            rop: ≈2^47.0, red: ≈2^47.0, δ: 1.011391, β: 61, d: 276, tag: euclidean

            >>> SIS.lattice(params.updated(norm=oo, length_bound=16), red_shape_model="lgsa")
            rop: ≈2^61.0, red: ≈2^59.9, sieve: ≈2^60.1, β: 95, η: 126, ζ: 0, d: 2486, prob: 1, ↻: 1, tag: infinity

            >>> SIS.lattice(params.updated(norm=oo, length_bound=16), red_shape_model="cn11")
            rop: ≈2^65.9, red: ≈2^64.9, sieve: ≈2^64.9, β: 113, η: 142, ζ: 0, d: 2486, prob: 1, ↻: 1, tag: infinity

        The success condition for euclidean norm bound is derived by determining the root hermite factor required for
        BKZ to produce the required output. For infinity norm bounds, the success conditions are derived using a
        probabilistic analysis. Vectors are assumed to be short as in [MATZOV22]_ P.18, or [Dilithium21]_ P.35.

        .. note :: When using euclidean norm bounds and the length bound is too small, this function returns
         β = d, and rop: inf

        """
        if params.norm == 2:
            tag = "euclidean"
        elif params.norm == oo:
            tag = "infinity"
        else:
            raise NotImplementedError(
                "SIS attack estimation currently only supports euclidean and infinity norms"
            )

        if tag == "infinity":
            red_shape_model = simulator_normalize(red_shape_model)

            f = partial(
                self.cost_zeta,
                params=params,
                red_shape_model=red_shape_model,
                red_cost_model=red_cost_model,
                log_level=log_level + 1,
            )

            if zeta is None:
                with local_minimum(0, params.m, log_level=log_level) as it:
                    for zeta in it:
                        it.update(
                            f(
                                zeta=zeta,
                                **kwds,
                            )
                        )
                # TODO: this should not be required
                cost = min(it.y, f(0, **kwds))
            else:
                cost = f(zeta=zeta)

        else:
            cost = self.cost_euclidean(
                params=params,
                red_cost_model=red_cost_model,
                log_level=log_level + 1,
            )

        cost["tag"] = tag
        cost["problem"] = params

        if tag == "euclidean":
            for k in ("sieve", "prob", "repetitions", "zeta"):
                try:
                    del cost[k]
                except KeyError:
                    pass

        return cost.sanity_check()

    __name__ = "lattice"


lattice = SISLattice()
