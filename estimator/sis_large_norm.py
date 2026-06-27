# -*- coding: utf-8 -*-
"""
Estimate cost of the attack of [C:DucEspPos23]_ on q-ary lattices (the small modulus regime).

See :ref:`SIS Lattice Attacks` for an introduction to what is available.

"""
from functools import partial
from math import log as mlog, sqrt as msqrt, exp as mexp, pi, erf

from scipy.optimize import brentq
from sage.all import oo, sqrt, RR, floor, cached_function

from .reduction import cost as costf
from .util import local_minimum, log2
from .cost import Cost
from .sis_parameters import SISParameters
from .simulator import normalize as simulator_normalize
from .prob import amplify as prob_amplify
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import max_beta

# tolerance, in length units, for identifying q-vectors when reading a simulated basis profile
profile_precision = 1e-6


def log2_lift_proportion_exact(n_q, sq_radius, q):
    """
    Exact base-2 logarithm of the proportion of ``Cube_{n_q}(q)`` lying inside a Euclidean ball.

    This counts the integer points directly, by convolving truncated theta series of ``ZZ`` as in
    [C:DucEspPos23]_ Sect. 3.3. It is a slow reference, of cost ``O(n_q · sq_radius · q)``, against
    which :func:`log2_lift_proportion` is validated; the latter should be preferred in practice.

    Both parities of ``q`` are handled: the centred representatives are ``{-(q-1)/2, ..., (q-1)/2}``
    for odd ``q`` and the non-symmetric ``{-(q-2)/2, ..., q/2}`` for even ``q`` ([C:DucEspPos23]_
    Def. 16, the latter being the case its footnote 3 sets aside).

    :param n_q: number of (uniform mod q) lifted coordinates.
    :param sq_radius: squared radius of the ball.
    :param q: the modulus.

    EXAMPLES::

        >>> from estimator.sis_large_norm import log2_lift_proportion_exact
        >>> log2_lift_proportion_exact(60, 40000, 127)  # doctest: +ELLIPSIS
        -19.99...

    """
    import numpy as np

    if sq_radius < 0:
        return RR(-oo)
    # the longest lift has every coordinate at the cube's outer edge, of magnitude q // 2
    max_sq = n_q * (q // 2) ** 2
    if sq_radius >= max_sq:
        # the whole cube lies inside the ball, every lift is short enough
        return RR(0)

    # squared-length distribution of one centred coordinate uniform mod q, truncated beyond X^D. The
    # range gives the q centred representatives for both parities: symmetric for odd q, and the
    # non-symmetric {-(q-2)/2, ..., q/2} for even q ([C:DucEspPos23]_ Def. 16).
    D = int(floor(sq_radius))
    dist = np.zeros(D + 1)
    for u in range(-((q - 1) // 2), q // 2 + 1):
        s = u * u
        if s <= D:
            dist[s] += 1.0 / q
    exps = np.nonzero(dist)[0]
    coeffs = dist[exps]

    # the n_q-fold convolution is the squared-length distribution of a lift; sum its mass in the ball
    p = np.zeros(D + 1)
    p[0] = 1.0
    for _ in range(n_q):
        nxt = coeffs[0] * p
        for j, c in zip(exps[1:], coeffs[1:]):
            nxt[int(j):] += c * p[: -int(j)]
        p = nxt
    total = float(p.sum())
    return RR(log2(total)) if total > 0 else RR(-oo)


def log2_lift_proportion(n_q, sq_radius, q):
    """
    Base-2 logarithm of the proportion of ``Cube_{n_q}(q)`` lying inside a Euclidean ball.

    ``Cube_{n_q}(q)`` is a centred set of representatives of ``(ZZ/qZZ)^{n_q}`` and the proportion of
    it inside the ball of squared radius ``sq_radius`` is the success probability of a single lift
    [C:DucEspPos23]_ Eq. 2. Rather than counting integer points directly, as in
    :func:`log2_lift_proportion_exact`, we apply the saddle-point (Bahadur-Rao) sharp large-deviation
    approximation [BahRao60]_ to the lower tail of the squared lift length. A single lifted coordinate
    is modelled by a random variable ``U`` continuous and uniform on ``[-q/2, q/2]`` -- giving each of
    the ``q`` integer representatives a unit-width tile, so the support has measure ``q`` and the same
    second moment up to ``O(1/q^2)`` for both parities of ``q``. The squared length ``U^2`` then has
    the closed-form cumulant generating function ``K(t) = log E[e^{t U^2}]``, so the estimate costs
    ``O(1)`` rather than the ``O(q)`` of the direct count and is accurate in the relevant (lower) tail
    (see the example below). The combination is, to our knowledge, not used elsewhere for this count,
    so the agreement with :func:`log2_lift_proportion_exact` is the warrant for correctness.

    :param n_q: number of (uniform mod q) lifted coordinates.
    :param sq_radius: squared radius of the ball.
    :param q: the modulus.

    EXAMPLES::

        >>> from estimator.sis_large_norm import log2_lift_proportion, log2_lift_proportion_exact
        >>> approx = log2_lift_proportion(60, 40000, 127)
        >>> exact = log2_lift_proportion_exact(60, 40000, 127)
        >>> bool(abs(approx - exact) < 0.1)
        True

    """
    a = float(sq_radius)
    h = q / 2.0  # a lifted coordinate is modelled as uniform on [-q/2, q/2], with E[U^2] = h^2 / 3
    if a >= n_q * h ** 2 / 3:
        # past the mean cube length the saddle-point's lower-tail regime ends; the proportion is ≥ 1/2
        # there and we round it up to 1, which is conservative and does not move the optimum, since the
        # number of lifts then saturates the per-attempt success probability (validated, see tests)
        return RR(0)
    if a < 1:
        # only the zero vector of the cube lies in the ball
        return RR(-n_q * log2(q))

    # below this value of the saddle parameter s = -t the closed-form integrals below degenerate to a
    # 0/0 form and lose accuracy to cancellation, so we use their s → 0 limit; this is a numerical
    # cutoff internal to the saddle point, unrelated to the q-vector tolerance ``profile_precision``
    small_s_cutoff = 1e-8

    def derivatives(t):
        # K'(t), K''(t) of K(t) = log E[e^{t U^2}] via I_k = ∫_0^h u^k e^{t u^2} du for one coordinate
        s = -t
        if s < small_s_cutoff:
            I0, I2, I4 = h, h ** 3 / 3, h ** 5 / 5
        else:
            decay = mexp(-s * h ** 2)
            I0 = msqrt(pi / s) / 2 * erf(h * msqrt(s))
            I2 = (I0 - h * decay) / (2 * s)
            I4 = (3 * I2 - h ** 3 * decay) / (2 * s)
        return I0, I2 / I0, I4 / I0 - (I2 / I0) ** 2

    # saddle point t < 0: the tilted mean n_q ⋅ K'(t) matches the squared radius
    t = brentq(lambda x: n_q * derivatives(x)[1] - a, -(2 * n_q / a + 50), 0.0)
    I0, _, var_t = derivatives(t)
    if var_t <= 0:
        return RR(0)
    ln_p = n_q * (mlog(I0) - mlog(h)) - t * a - mlog(-t) - mlog(2 * pi * n_q * var_t) / 2
    return RR(min(0.0, ln_p) / mlog(2))


class SISLargeNorm:
    """
    Estimate cost of the attack of [C:DucEspPos23]_ on q-ary lattices (the small modulus regime).
    """

    @staticmethod
    @cached_function
    def cost(
        beta: int,
        params: SISParameters,
        success_probability: float = 0.99,
        red_shape_model="ZGSA",
        red_cost_model=red_cost_model_default,
        exact: bool = False,
        log_level=None,
        **kwds,
    ):
        """
        Cost of the attack of [C:DucEspPos23]_ on q-ary lattices for a fixed block size β.

        The attack reduces a basis of the kernel lattice ``Λ_q^⊥(A)``, which exhibits a "Z-shape":
        a head of ``n_q`` vectors of length ``q``, a sloped middle and a flat tail. Sieving in the
        projected sublattice past the q-vectors yields many short vectors which are then lifted over
        the q-vectors; a lift is a solution when its (uniform mod q) lifted entries are short enough.
        We take the cheaper of lifting the terminal sieve database [C:DucEspPos23]_ Sect. 3.4 and
        lifting every vector of the final sieve iteration [C:DucEspPos23]_ Sect. 4.

        :param beta: block size β.
        :param params: SIS parameters.
        :param success_probability: targeted success probability.
        :param red_shape_model: how to model the shape of a reduced basis; must expose the q-vectors.
        :param red_cost_model: how to cost lattice reduction.
        :param exact: count the short lifts with the slow exact :func:`log2_lift_proportion_exact`
            instead of the fast saddle-point :func:`log2_lift_proportion` (default).

        """
        d = params.m

        if beta > d:
            return Cost(rop=oo)

        # shape of the BKZ-β reduced basis of the kernel lattice Λ_q^⊥(A)
        simulator = simulator_normalize(red_shape_model)
        r = simulator(d=d, n=d - params.n, q=params.q, beta=beta, xi=1, tau=False)

        # number of q-vectors n_q left at the head of the basis (Zone I of the Z-shape)
        if abs(sqrt(r[0]) - params.q) < profile_precision:  # q-vectors exist
            n_q = next((i for i, r_ in enumerate(r) if r_ < r[0]), len(r))
        else:
            # no q-vectors to lift over: this is the plain reduction attack of SIS.lattice
            return Cost(rop=oo)

        # the sieve runs in the projected sublattice Λ_{[ℓ:s]}, ℓ = n_q + 1 and s = min(ℓ+β, d+1)
        sieve_dim = min(beta, d - n_q)

        # the fast saddle-point count by default; the slow exact count when opted in
        lift_proportion = log2_lift_proportion_exact if exact else log2_lift_proportion

        # expected number of solutions from one sieve over the two lifting strategies; for each the
        # projected vectors have length ``proj_length`` and there are ``2^(log2_density ⋅ sieve_dim)``
        log2_solutions = None
        for proj_length, log2_density in (
            (sqrt(RR(4) / 3) * params.q, log2(RR(4) / 3) / 2),  # terminal sieve database, Sect. 3.4
            (sqrt(RR(2)) * params.q, log2(RR(3) / 2) / 2),  # on-the-fly lifting, Sect. 4
        ):
            if proj_length >= params.length_bound:
                # the projected vectors alone are already too long to lift to a solution
                continue
            # squared lifting radius: the lifted entries must absorb the remaining length, Eq. 2
            sq_radius = floor(RR(params.length_bound) ** 2 - proj_length ** 2)
            log2_p = log2_density * sieve_dim + lift_proportion(n_q, sq_radius, params.q)
            log2_solutions = log2_p if log2_solutions is None else max(log2_solutions, log2_p)

        if log2_solutions is None:
            return Cost(rop=oo)

        probability = RR(2) ** min(RR(0), log2_solutions)

        # cost of producing the short vectors and the underlying BKZ-β reduction
        rho, cost_red, N, sieve_dim = red_cost_model.short_vectors(beta, d, sieve_dim=sieve_dim)
        bkz_cost = costf(red_cost_model, beta, d)

        ret = Cost()
        ret["rop"] = cost_red
        ret["red"] = bkz_cost["rop"]
        ret["sieve"] = max(cost_red - bkz_cost["rop"], 1e-100)  # ensuring non-zero cost here
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
        # repeat the whole attack ~1/prob times, rerandomising the basis between attempts
        if probability and not RR(probability).is_NaN():
            ret = ret.repeat(prob_amplify(success_probability, probability))
        else:
            return Cost(rop=oo)

        return ret

    def __call__(
        self,
        params: SISParameters,
        red_shape_model="ZGSA",
        red_cost_model=red_cost_model_default,
        exact: bool = False,
        log_level=1,
        **kwds,
    ):
        """
        Estimate the cost of the attack of [C:DucEspPos23]_ on q-ary lattices.

        :param params: SIS parameters.
        :param red_shape_model: how to model the shape of a reduced basis; must expose the q-vectors.
        :param red_cost_model: how to cost lattice reduction.
        :param exact: count the short lifts with the slow exact :func:`log2_lift_proportion_exact`
            instead of the fast saddle-point :func:`log2_lift_proportion` (default). This is much
            slower but removes the approximation, for when one wants certainty and can wait.
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
            >>> SIS.large_norm(params)
            rop: ≈2^46.6, red: ≈2^46.6, sieve: ≈2^40.3, β: 59, η: 59, ℓ: 32, d: 300, prob: 1, ↻: 1, tag: large_norm

        Setting ``exact`` replaces the fast saddle-point count by the slow exact one; it is much
        slower but here confirms the default::

            >>> SIS.large_norm(params, exact=True)
            rop: ≈2^46.6, red: ≈2^46.6, sieve: ≈2^40.3, β: 59, η: 59, ℓ: 32, d: 300, prob: 1, ↻: 1, tag: large_norm

        The attack only targets the regime where the length bound exceeds the modulus and the
        euclidean norm is used; outside of it ``rop`` is infinite::

            >>> SIS.large_norm(params.updated(length_bound=200))
            rop: ≈2^inf, tag: large_norm

        """
        if params.norm != 2:
            raise NotImplementedError(
                "This attack is only available for the euclidean norm."
            )

        # the attack targets the regime where the length bound exceeds the modulus [C:DucEspPos23]_
        if params.length_bound <= params.q:
            cost = Cost(rop=oo)
            cost["tag"] = "large_norm"
            cost["problem"] = params
            return cost

        f = partial(
            SISLargeNorm.cost,
            params=params,
            red_shape_model=red_shape_model,
            red_cost_model=red_cost_model,
            exact=exact,
        )

        with local_minimum(40, min(params.m, max_beta), precision=2, log_level=log_level) as it:
            for beta in it:
                it.update(f(beta=beta))
            for beta in it.neighborhood:
                it.update(f(beta=beta))
            cost = it.y

        if cost is None:
            cost = Cost(rop=oo)

        Logging.log("large_norm", log_level, f"{cost!r}")

        cost["tag"] = "large_norm"
        cost["problem"] = params
        return cost.sanity_check()

    __name__ = "large_norm"


large_norm = SISLargeNorm()
