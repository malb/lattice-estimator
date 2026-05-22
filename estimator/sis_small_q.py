# -*- coding: utf-8 -*-
"""
Estimate cost of solving SIS with a small modulus.

See :ref:`SIS Lattice Attacks` for an introduction to what is available.

"""
from functools import partial
from math import log as mlog, pi

import numpy as np
from scipy.optimize import brentq
from sage.all import oo, sqrt, RR, floor, cached_function

from .reduction import cost as costf
from .util import local_minimum, log2
from .cost import Cost
from .sis_parameters import SISParameters
from .simulator import GSA, normalize as simulator_normalize
from .prob import amplify as prob_amplify
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import red_shape_model as red_shape_model_default
from .conf import max_beta


@cached_function
def _centred_squares(q):
    """
    Squares of a centred set of representatives of ``ZZ/qZZ``.

    These are the squared lengths a single (uniform mod q) lifted coordinate can take, used to form
    the cumulant generating function to the squared length of a lift.

    :param q: the modulus.

    """
    u = np.arange(-((q - 1) // 2), q // 2 + 1, dtype=float)
    return u * u


def log2_lift_proportion(n_q, sq_radius, q):
    """
    Base-2 logarithm of the proportion of ``Cube_{n_q}(q)`` lying inside a Euclidean ball.

    ``Cube_{n_q}(q)`` is a centred set of representatives of ``(ZZ/qZZ)^{n_q}`` and the proportion of
    it inside the ball of squared radius ``sq_radius`` is the success probability of a single lift
    [C:DucEspPos23]_ Eq. 2. We evaluate it with a saddle-point (Bahadur-Rao) approximation of the
    theta-series count of [C:DucEspPos23]_ Sect. 3.3, in closed form and accurate in the relevant
    tail. Above the mean cube length the lifting is not the bottleneck and the proportion is one.

    :param n_q: number of (uniform mod q) lifted coordinates.
    :param sq_radius: squared radius of the ball.
    :param q: the modulus.

    """
    squares = _centred_squares(q)
    if sq_radius >= n_q * float(squares.mean()):
        return RR(0)
    if sq_radius < 1:
        # only the zero vector of the cube lies in the ball
        return RR(-n_q * mlog(q) / mlog(2))

    def derivatives(t):
        # φ(t) and the derivatives K'(t), K''(t) of K(t) = log φ(t) for one coordinate's squared length
        w = np.exp(t * squares)
        s0, s1, s2 = float(w.sum()), float((squares * w).sum()), float((squares * squares * w).sum())
        return s0, s1 / s0, s2 / s0 - (s1 / s0) ** 2

    # saddle point t < 0: the tilted mean n_q ⋅ K'(t) matches the squared radius
    t = brentq(lambda t: n_q * derivatives(t)[1] - sq_radius, -50.0, 0.0)
    s0, _, var_t = derivatives(t)
    log_phi = mlog(s0) - mlog(q)
    ln_p = n_q * log_phi - t * sq_radius - mlog(-t) - mlog(2 * pi * n_q * var_t) / 2
    return RR(min(0.0, ln_p) / mlog(2))


class SISSmallQ:
    """
    Estimate cost of solving SIS with a small modulus via lattice reduction.
    """

    @staticmethod
    @cached_function
    def cost(
        beta: int,
        params: SISParameters,
        on_the_fly: bool = None,
        inhom: str = None,
        success_probability: float = 0.99,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=None,
        **kwds,
    ):
        """
        Cost of the small-q attack [C:DucEspPos23] for a fixed block size β.

        The attack reduces a basis of the kernel lattice ``Λ_q^⊥(A)``, which exhibits a "Z-shape":
        a head of ``n_q`` vectors of length ``q``, a sloped middle and a flat tail. Sieving in the
        projected sublattice past the q-vectors yields many short vectors which are then lifted over
        the q-vectors; a lift is a solution when its (uniform mod q) lifted entries are short enough.

        :param beta: block size β.
        :param params: SIS parameters.
        :param on_the_fly: lift every vector seen in the last sieve iteration [C:DucEspPos23]_ Sect. 4.
            If ``None``, return the cheaper of the two variants.
        :param inhom: model the reduction to the inhomogeneous problem [C:DucEspPos23]_ Sect. 3.5;
            one of ``None`` (homogeneous SIS*), ``"specific"`` (loss ``q/2``) or ``"generic"``
            (loss ``m(q-1)``).
        :param success_probability: targeted success probability.
        :param red_shape_model: how to model the shape of a reduced basis; must expose the q-vectors.
        :param red_cost_model: how to cost lattice reduction.

        """
        if on_the_fly is None:
            # the two lifting strategies trade more vectors for a smaller radius; consider the cheaper
            return min(
                SISSmallQ.cost(beta, params, otf, inhom, success_probability,
                               red_shape_model, red_cost_model, log_level)
                for otf in (True, False)
            )

        # the inhomogeneous reduction operates one rank higher
        d = params.m + 1 if inhom else params.m

        if beta > d:
            return Cost(rop=oo)

        # length and (log) number of vectors expected from sieving the projected sublattice
        if on_the_fly:
            # lift every vector visited in the final sieve iteration 
            proj_length = sqrt(RR(2)) * params.q
            log2_density = log2(RR(3) / 2) / 2
        else:
            # lift the terminal sieve database only
            proj_length = sqrt(RR(4) / 3) * params.q
            log2_density = log2(RR(4) / 3) / 2

        if proj_length >= params.length_bound:
            # the projected vectors alone are already too long to lift to a solution
            return Cost(rop=oo)

        # shape of the BKZ-β reduced basis of the kernel lattice Λ_q^⊥(A)
        simulator = simulator_normalize(red_shape_model)
        if simulator is GSA:
            # the GSA has no q-vectors; the attack requires the Z-shape, so use the ZGSA instead
            simulator = simulator_normalize("ZGSA")
        r = simulator(d=d, n=d - params.n, q=params.q, beta=beta, xi=1, tau=False)

        # number of q-vectors n_q left at the head of the basis (Zone I of the Z-shape)
        if abs(sqrt(r[0]) - params.q) < 1e-6:  # q-vectors exist
            n_q = next((i for i, r_ in enumerate(r) if r_ < r[0]), len(r))
        else:
            # no q-vectors to lift over: this is the plain reduction attack of SIS.lattice
            return Cost(rop=oo)

        # the sieve runs in the projected sublattice Λ_{[ℓ:s]}, ℓ = n_q + 1 and s = min(ℓ+β, d+1)
        sieve_dim = min(beta, d - n_q)

        # squared lifting radius: the lifted entries must absorb the remaining length 
        sq_radius = floor(RR(params.length_bound) ** 2 - proj_length ** 2)
        log2_p = log2_lift_proportion(n_q, sq_radius, params.q)

        # expected number of solutions from one sieve, capped at one 
        log2_solutions = log2_density * sieve_dim + log2_p

        # the inhomogeneous reductions lose a factor in success probability 
        if inhom == "specific":
            log2_solutions += log2(RR(2) / params.q)
        elif inhom == "generic":
            log2_solutions += -log2((params.m + 1) * (params.q - 1))
        elif inhom is not None:
            raise ValueError(f"Unknown inhomogeneous reduction: {inhom}")

        probability = RR(2) ** min(RR(0), log2_solutions)

        # cost of producing the short vectors and the underlying BKZ-β reduction
        rho, cost_red, N, sieve_dim = red_cost_model.short_vectors(beta, d, sieve_dim=sieve_dim)
        bkz_cost = costf(red_cost_model, beta, d)

        ret = Cost()
        ret["rop"] = cost_red
        ret["red"] = bkz_cost["rop"]
        ret["sieve"] = max(cost_red - bkz_cost["rop"], 1e-100)  # check non-zero cost
        ret["beta"] = beta
        ret["eta"] = sieve_dim
        ret["ell"] = n_q + 1
        ret["d"] = d
        ret["prob"] = probability
        ret.register_impermanent(
            rop=True,
            red=True,
            sieve=True,
            eta=False,
            ell=False,
            prob=False,
        )
        # repeat the whole attack ~1/prob times, rerandomizing the basis between attempts
        if probability and not RR(probability).is_NaN():
            ret = ret.repeat(prob_amplify(success_probability, probability))
        else:
            return Cost(rop=oo)

        return ret

    def __call__(
        self,
        params: SISParameters,
        on_the_fly: bool = None,
        inhom: str = None,
        red_shape_model=red_shape_model_default,
        red_cost_model=red_cost_model_default,
        log_level=1,
        **kwds,
    ):
        """
        Estimate the cost of solving SIS with a small modulus via reduction.

        :param params: SIS parameters.
        :param on_the_fly: lift every vector seen in the last sieve iteration [C:DucEspPos23]_ Sect. 4.
            If ``None``, report the cheaper of the two variants.
        :param inhom: model the reduction to the inhomogeneous problem [C:DucEspPos23]_ Sect. 3.5;
            one of ``None`` (homogeneous SIS*), ``"specific"`` (loss ``q/2``) or ``"generic"``
            (loss ``m(q-1)``).
        :param red_shape_model: how to model the shape of a reduced basis; must expose the q-vectors.
        :param red_cost_model: how to cost lattice reduction.
        :return: A cost dictionary.

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction.
        - ``sieve``: Number of word operations in sieving.
        - ``β``: BKZ block size.
        - ``η``: Dimension of the sieve in the projected sublattice.
        - ``ℓ``: Index of the projected sublattice, one more than the number of q-vectors.
        - ``d``: Lattice dimension.
        - ``prob``: Probability of success of one attempt.
        - ``repeat``: How often to repeat the attack.

        EXAMPLES::

            >>> from estimator import *
            >>> params = SIS.Parameters(n=150, q=257, length_bound=420, m=300, norm=2)
            >>> SIS.small_q(params)
            rop: ≈2^46.6, red: ≈2^46.6, sieve: ≈2^40.3, β: 59, η: 59, ℓ: 32, d: 300, prob: 1, ↻: 1, tag: small_q

        The reduction to the inhomogeneous problem (e.g. signature forgery) loses a factor in the
        success probability and operates one rank higher [C:DucEspPos23]_ Sect. 3.5::

            >>> SIS.small_q(params, inhom="specific")
            rop: ≈2^48.3, red: ≈2^48.2, sieve: ≈2^42.1, β: 65, η: 65, ℓ: 27, d: 301, prob: 1, ↻: 1, tag: small_q_isis

        The attack only targets the regime where the length bound exceeds the modulus and the
        euclidean norm is used; outside of it ``rop`` is infinite::

            >>> SIS.small_q(params.updated(length_bound=200))
            rop: ≈2^inf, tag: small_q

        """
        if params.norm != 2:
            raise NotImplementedError(
                "The small-q attack is only available for the euclidean norm."
            )

        # the attack targets the regime ν > q 
        if params.length_bound <= params.q:
            cost = Cost(rop=oo)
            cost["tag"] = "small_q"
            cost["problem"] = params
            return cost

        f = partial(
            SISSmallQ.cost,
            params=params,
            on_the_fly=on_the_fly,
            inhom=inhom,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
        )

        with local_minimum(40, min(params.m, max_beta), precision=2, log_level=log_level) as it:
            for beta in it:
                it.update(f(beta=beta))
            for beta in it.neighborhood:
                it.update(f(beta=beta))
            cost = it.y

        if cost is None:
            cost = Cost(rop=oo)

        Logging.log("sis_small_q", log_level, f"{cost!r}")

        cost["tag"] = "small_q" if inhom is None else "small_q_isis"
        cost["problem"] = params
        return cost.sanity_check()

    __name__ = "small_q"


small_q = SISSmallQ()
