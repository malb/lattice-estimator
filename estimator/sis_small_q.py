# -*- coding: utf-8 -*-
"""
Estimate cost of small modulus (I)SIS attack [DEP23].

See :ref:`SIS Small-q Attacks` for an introduction what is available

"""
from sage.all import oo, RR, ceil
from .cost import Cost
from .io import Logging

# for speed in hot paths and fast vectorized convolutions
from math import log as mlog, sqrt as msqrt, floor as mfloor, exp as mexp
from numpy import zeros as np_zeros

LOG_INFINITY = 9999

# boxed theta: counting integer points in ball-cube intersection
class BoxedTheta:

    def __init__(self, R, q):
        """
        :param R: the radius of the ball
        :param q: the side length of the box (modulus)
        """
        self.ell = int(mfloor(R**2 + 1))

        # theta0 has all probability on length 0
        Theta0 = np_zeros(self.ell)
        Theta0[0] = 1.0

        # theta1 gives probability distribution over squared lengths for uniform mod q
        Theta1 = np_zeros(self.ell)
        for i in range(q):
            x = min(i, q - i)
            x2 = x * x
            if x2 < self.ell:
                Theta1[x2] += 1.0 / q

        self.cache = {0: Theta0, 1: Theta1}

    def _convol(self, v, w):
        """
        convolve two theta series by numpy vectorized operations
        """
        x = v[0] * w.copy()
        for i in range(1, self.ell):
            if v[i] == 0:
                continue
            x[i:] += float(v[i]) * w[:-i]
        return x

    def __call__(self, n):
        """
        compute the probability distribution over squared lengths for n dimensions.

        :param n: the number of uniform mod q coordinates
        :returns: numpy array with probability of each squared length < self.ell
        """
        if n in self.cache:
            return self.cache[n]

        # recursively build up
        cached = self(n - 1)
        x = self._convol(self.cache[1], cached)
        self.cache[n] = x
        return x


# cache for boxed theta instances
_bt_cache = {}


def log_intersection_proportion(n, R, q):
    """
    compute the natural log of the probability that a uniform vector mod q
    (centered around 0) falls into an l2 ball of radius R.

    :returns: natural log of the probability
    """
    if R <= 0:
        return float(-LOG_INFINITY)

    key = (float(R), int(q))
    if key in _bt_cache:
        bt = _bt_cache[key]
    else:
        bt = BoxedTheta(R, q)
        _bt_cache[key] = bt

    proportion = bt(n)
    total = float(proportion.sum())
    if total <= 0:
        return float(-LOG_INFINITY)
    return mlog(total)



# BKZ shape modeling for q-ary lattices (z-shape)
def delta_0f(k):
    """
    root Hermite factor for BKZ-k; small values are experimentally determined, otherwise from [Che13].

    :param k: BKZ block size
    :returns: root Hermite factor delta_0
    """
    from math import pi, e

    small = (
        (2, 1.02190),
        (5, 1.01862),
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    )

    k = float(k)
    if k <= 2:
        return 1.0219
    elif k < 40:
        for i in range(1, len(small)):
            if small[i][0] > k:
                return small[i - 1][1]
    elif k == 40:
        return small[-1][1]
    else:
        return (k / (2 * pi * e) * (pi * k) ** (1.0 / k)) ** (1.0 / (2 * (k - 1.0)))


def svp_classical(b):
    """
    log_2 of best known classical cost of SVP in dimension b; use the exponent 0.292 from [BDGL16].

    :param b: SVP dimension
    :returns: log_2 of classical SVP cost
    """
    from math import log as mlog, sqrt as msqrt
    return b * mlog(msqrt(3.0 / 2)) / mlog(2)  # 0.292 * b


def construct_bkz_shape(q, nq, n1, b):
    """
    simulate the (log) shape of a basis after BKZ-b reduction on a q-ary lattice

    the initial basis has nq q-vectors and n1 unit vectors. the attack uses the Z-shape: q-vectors remain at the start, followed by a GSA slope, then unit vectors.

    :param q: the modulus
    :param nq: initial number of q-vectors
    :param n1: number of unit vectors (in the flat tail)
    :param b: BKZ block size
    :returns: tuple (a, a+B, L) where:
        - a: number of remaining q-vectors
        - a+B: end index of the sloped region
        - L: the log-profile of the basis
    """
    d = nq + n1
    lq = mlog(q)
    glv = nq * lq  # goal log volume

    slope = -2 * mlog(delta_0f(b))

    if b == 0:
        L = [lq] * nq + [0.0] * n1
        return (nq, nq, L)

    B = int(mfloor(lq / (-slope)))  # number of vectors in the sloped region
    L = [lq] * nq + [lq + i * slope for i in range(1, B + 1)] + [0.0] * n1

    x = 0
    lv = sum(L[:d])

    # slide window right until volume matches goal
    while lv > glv and x + d < len(L):
        lv -= L[x]
        lv += L[x + d]
        x += 1

    if x > B:
        x = B

    L = L[x:x + d]
    a = max(0, nq - x)  # length of [q, ..., q] sequence
    B = min(B, d - a)   # length of GSA sequence

    # small shift to balance/equilibrate volume
    diff = glv - lv
    if B > 0 and abs(diff) >= 1e-10:
        for i in range(a, a + B):
            L[i] += diff / B

    return (a, a + B, L)


def construct_bkz_shape_randomised(q, nq, n1, b):
    """
    simulate the shape after randomisation (killing q-vectors) and BKZ-b reduction

    :returns: tuple (a, a+B, L) as in construct_bkz_shape
    """
    glv = nq * mlog(q)
    d = nq + n1
    L = []

    slope = -2 * mlog(delta_0f(b))
    li = 0.0
    lv = 0.0

    for _ in range(d):
        li -= slope
        lv += li
        if lv > glv:
            break
        L = [li] + L

    B = len(L)
    L += [0.0] * (d - B)
    a = 0

    lv = sum(L)
    diff = lv - glv
    if B > 0:
        for i in range(a, a + B):
            L[i] -= diff / B

    return (a, a + B, L)


def bkz_first_length(q, nq, n1, b):
    """
    estimate the first vector length after randomisation and BKZ-b reduction.

    :returns: estimated first vector length
    """
    (_, _, L) = construct_bkz_shape_randomised(q, nq, n1, b)
    return mexp(L[0])



# attack cost estimation
def success_prob(q, q_block_len, nu, proj_len):
    """
    compute the log probability that a single lift succeeds; a lift succeeds if the lifted vector has length <= nu

    :param q_block_len: number of q-vectors to lift over
    :param nu: length bound
    :param proj_len: length of the projected vector
    :returns: natural log of success probability
    """
    if q_block_len <= 0:
        return float(-LOG_INFINITY)
    if proj_len > nu:
        return float(-LOG_INFINITY)

    R = msqrt(nu**2 - proj_len**2)
    return log_intersection_proportion(q_block_len, R, q)


def sis_l2_cost_randomised(q, w, h, nu, b, cost_svp=svp_classical, verbose=False):
    """
    cost of standard lattice reduction attack (randomizing to remove q-vectors)

    :param w: lattice rank (m)
    :param h: number of equations (n), so volume = q^h
    :param cost_svp: SVP cost function
    :returns: log_2 cost
    """
    first_len = bkz_first_length(q, h, w - h, b)
    if first_len > nu:
        return float(LOG_INFINITY)
    return cost_svp(b)


def sis_l2_cost_small_q(q, w, h, nu, b, cost_svp=svp_classical, verbose=False,
                        sieve="bdgl", otf_lift=True, inhom=None):
    """
    cost of the small-q attack, which exploits Z-shape

    :param sieve: sieve type ("svp_only", "nv", "bdgl")
    :param otf_lift: whether to use on-the-fly lifting
    :param inhom: None for SIS*, "worst" for generic ISIS, "specific" for optimized ISIS
    :returns: log_2 cost
    """
    if sieve not in {"svp_only", "nv", "bdgl"}:
        raise ValueError(f"Unknown sieve type: {sieve}")
    if inhom is not None and inhom not in {"worst", "specific"}:
        raise ValueError(f"Unknown inhom type: {inhom}")

    # for ISIS, add one dimension for the reduction
    if inhom in {"worst", "specific"}:
        w += 1

    (q_block_len, _, _) = construct_bkz_shape(q, h, w - h, b)
    first_len = q

    if q_block_len == 0 or sieve == "svp_only":
        # BKZ already removed q-vector structure
        return sis_l2_cost_randomised(q, w, h, nu, b, cost_svp=cost_svp, verbose=verbose)

    if not otf_lift:
        # saturated ball of radius sqrt(4/3)
        proj_len = first_len * msqrt(4.0 / 3.0)
        log_projs = b * mlog(4.0 / 3.0) / 2.0
    else:
        if sieve == "nv":
            # nguyen-vidick sieve: more vectors, longer
            proj_len = first_len * msqrt(4.0 / 3.0) * msqrt(2.0)
            log_projs = b * mlog(4.0 / 3.0)
        elif sieve == "bdgl":
            # [BDGL16] sieve: fewer vectors, shorter
            proj_len = first_len * msqrt(2.0)
            log_projs = b * mlog(3.0 / 2.0) / 2.0

    if proj_len > nu:
        return float(LOG_INFINITY)

    log_success = success_prob(q, q_block_len, nu, proj_len)

    if inhom is None:
        # solve SIS*
        log_p = min(0, log_success + log_projs)
    elif inhom == "worst":
        # solve ISIS with generic reduction (probability loss ~mq)
        log_p = min(0, log_success + log_projs - mlog(w * (q - 1)))
    elif inhom == "specific":
        # solve ISIS with optimized reduction (probability loss q/2)
        log_p = min(0, log_success + log_projs + mlog(2.0 / q))

    cost = cost_svp(b) - log_p / mlog(2.0)

    if verbose:
        Logging.log("small_q", 1, f"BKZ blocksize = {b}")
        Logging.log("small_q", 1, f"lift dim = {q_block_len}")
        Logging.log("small_q", 1, f"R_qlift = {round(float(msqrt(nu**2 - proj_len**2)), 5)}")
        Logging.log("small_q", 1, f"log2 Success = {round(float(log_success / mlog(2.0)), 5)}")
        Logging.log("small_q", 1, f"log2_#targets = {round(float(log_projs / mlog(2)), 5)}")
        Logging.log("small_q", 1, f"log2_p = {round(float(log_p / mlog(2)), 5)}")
        Logging.log("small_q", 1, f"cost = {round(float(cost), 5)}")

    return cost


def optimize_small_q_attack(q, w, h, nu, cost_svp=svp_classical,
                            sieve="bdgl", otf_lift=True, inhom=None,
                            log_level=5):
    """
    find optimal BKZ block size for the small-q attack

    :param log_level: logging verbosity
    :returns: tuple (best_beta, best_cost)
    """
    best_cost = float(LOG_INFINITY)
    best_b = w

    # iterate from w down, with early abort
    for b in range(w, 2, -1):
        cost = sis_l2_cost_small_q(q, w, h, nu, b, cost_svp=cost_svp,
                                   sieve=sieve, otf_lift=otf_lift, inhom=inhom)
        if cost <= best_cost:
            best_cost = cost
            best_b = b
            Logging.log("small_q", log_level, f"b = {b}, cost = {cost:.3f}")

        if cost >= best_cost + 10:
            # found global minimum
            break

    return (best_b, best_cost)


# estimator interface
class SISSmallQ:

    @staticmethod
    def _is_applicable(params):
        # check if the small-q attack is likely beneficial; the attack is most effective when nu > q, but can also help when nu ~ q
        # the attack is relevant when nu is comparable to or larger than q
        return params.length_bound >= params.q / 2

    def __call__(
        self,
        params,
        sieve="bdgl",
        otf_lift=True,
        inhom=None,
        log_level=5,
    ):
        """
        estimate the cost of the small-q attack on SIS/ISIS

        :param params: SIS parameters
        :returns: Cost dictionary

        EXAMPLES::

            >>> from estimator import SIS
            >>> params = SIS.Parameters(n=512, q=257, length_bound=801, m=1024, norm=2)
            >>> SIS.small_q(params, otf_lift=False)  # doctest: +SKIP
            rop: ≈2^95.0, β: 325, d: 1024, prob: ≈2^-65.9, ↻: ≈2^65.9, tag: small_q_sis

            >>> SIS.small_q(params, otf_lift=True)  # doctest: +SKIP
            rop: ≈2^90.0, β: 308, d: 1024, prob: ≈2^-60.9, ↻: ≈2^60.9, tag: small_q_sis

            >>> SIS.small_q(params, inhom="specific", otf_lift=True)  # doctest: +SKIP
            rop: ≈2^92.0, β: 314, d: 1024, prob: ≈2^-6.0, ↻: ≈2^6.0, tag: small_q_isis

        """
        q = int(params.q)
        m = int(params.m)  # lattice rank w
        n = int(params.n)  # number of equations h
        nu = float(params.length_bound)

        Logging.log("small_q", log_level, f"small-q attack on q={q}, m={m}, n={n}, nu={nu}")

        # optimize block size
        best_b, best_cost = optimize_small_q_attack(
            q, m, n, nu,
            cost_svp=svp_classical,
            sieve=sieve,
            otf_lift=otf_lift,
            inhom=inhom,
            log_level=log_level + 1,
        )

        if best_cost >= LOG_INFINITY:
            return Cost(rop=oo)

        # compute success probability for reporting
        if inhom in {"worst", "specific"}:
            w = m + 1
        else:
            w = m

        (q_block_len, _, _) = construct_bkz_shape(q, n, w - n, best_b)

        if not otf_lift:
            proj_len = q * msqrt(4.0 / 3.0)
            log_projs = best_b * mlog(4.0 / 3.0) / 2.0
        else:
            if sieve == "bdgl":
                proj_len = q * msqrt(2.0)
                log_projs = best_b * mlog(3.0 / 2.0) / 2.0
            else:
                proj_len = q * msqrt(4.0 / 3.0) * msqrt(2.0)
                log_projs = best_b * mlog(4.0 / 3.0)

        log_success = success_prob(q, q_block_len, nu, proj_len)

        if inhom is None:
            log_p = min(0, log_success + log_projs)
        elif inhom == "worst":
            log_p = min(0, log_success + log_projs - mlog(w * (q - 1)))
        else:  # inhom == "specific"
            log_p = min(0, log_success + log_projs + mlog(2.0 / q))

        prob = RR(2) ** (log_p / mlog(2))
        repetitions = 1 / prob if prob > 0 else oo

        # set tag based on problem type
        if inhom is None:
            tag = "small_q_sis"
        else:
            tag = "small_q_isis"

        ret = Cost()
        ret["rop"] = RR(2) ** best_cost
        ret["beta"] = best_b
        ret["d"] = m
        ret["prob"] = prob
        ret["tag"] = tag
        ret["problem"] = params

        Cost.register_impermanent(beta=False, d=False, prob=False, tag=False)

        # apply repetition
        if prob > 0 and repetitions < oo:
            ret = ret.repeat(ceil(repetitions))

        return ret

    __name__ = "small_q"


small_q = SISSmallQ()