# -*- coding: utf-8 -*-
"""
Simulate lattice reduction on the rows of::

    ⌜ ξI  A  0 ⌝
    ǀ  0 qI  0 |
    ⌞ 0   c  τ ⌟

where

- ξI ∈ ZZ^{n × n},
- A ∈ ZZ_q^{n × m},
- qI ∈ ZZ^{m × m},
- τ ∈ ZZ and
- d = m + n + 1.

The last row is optional.
"""

from sage.all import RR, log, line, cached_function, pi, exp
from functools import partial


def qary_simulator(f, d, n, q, beta, xi=1, tau=1, dual=False, ignore_qary=False):
    """
    Reduced lattice shape calling ``f``.

    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n-1`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: Kannan factor τ.
    :param dual: perform reduction on the dual.
    :param ignore_qary: Ignore the special q-ary structure (forget q vectors)

    """

    if not tau:
        r = [q**2] * (d - n) + [xi**2] * n
    else:
        r = [q**2] * (d - n - 1) + [xi**2] * n + [tau**2]

    if ignore_qary:
        r = GSA(d, n, q, 2, xi=xi, tau=tau)

    if dual:
        # 1. reverse and reflect the basis (go to dual)
        r = [1 / r_ for r_ in reversed(r)]
        # 2. simulate reduction on the dual basis
        r = f(r, beta)
        # 3. reflect and reverse the basis (go back to primal)
        r = [1 / r_ for r_ in reversed(r)]
        return r
    else:
        return f(r, beta)


def CN11(d, n, q, beta, xi=1, tau=1, dual=False, ignore_qary=False):
    """
    Reduced lattice shape using simulator from [AC:CheNgu11]_

    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n-1`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: Kannan factor τ.
    :param dual: perform reduction on the dual.
    :param ignore_qary: Ignore the special q-ary structure (forget q vectors)
    :returns: squared Gram-Schmidt norms

    """

    from fpylll import BKZ
    from fpylll.tools.bkz_simulator import simulate

    def f(r, beta):
        return simulate(r, BKZ.EasyParam(beta))[0]

    return qary_simulator(f=f, d=d, n=n, q=q, beta=beta, xi=xi, tau=tau, dual=dual, ignore_qary=ignore_qary)


def GSA(d, n, q, beta, xi=1, tau=1, dual=False):
    """
    Reduced lattice shape following the Geometric Series Assumption [Schnorr03]_

    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n-1`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: Kannan factor τ.
    :param dual: ignored, since GSA is self-dual: applying the GSA to the dual is equivalent to
           applying it to the primal.
    :returns: squared Gram-Schmidt norms

    """
    from .reduction import delta as deltaf

    if not tau:
        log_vol = RR(log(q, 2) * (d - n) + log(xi, 2) * n)
    else:
        log_vol = RR(log(q, 2) * (d - n - 1) + log(xi, 2) * n + log(tau, 2))

    delta = deltaf(beta)
    r_log = [(d - 1 - 2 * i) * RR(log(delta, 2)) + log_vol / d for i in range(d)]
    r = [2 ** (2 * r_) for r_ in r_log]
    return r


def ZGSA(d, n, q, beta, xi=1, tau=1, dual=False):
    from math import lgamma
    from .util import gh_constant, small_slope_t8
    """
    Reduced lattice Z-shape following the Geometric Series Assumption as specified in
    NTRU fatrigue [DucWoe21]_
    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param dual: ignored, since GSA is self-dual: applying the GSA to the dual is equivalent to
           applying it to the primal.
    :returns: Squared Gram-Schmidt norms

    EXAMPLES:
        >>> from estimator.simulator import GSA, ZGSA, CN11
        >>> n = 128
        >>> d = 213
        >>> q = 2048
        >>> beta = 40
        >>> xi = 1
        >>> tau = 1
        >>> zgsa_profile = ZGSA(d, n, q, beta, xi, tau)
        >>> len(zgsa_profile)
        214

    Setting tau to False indicates a homogeneous instance.

        >>> tau = False
        >>> zgsa_profile = ZGSA(d, n, q, beta, xi, tau)
        >>> len(zgsa_profile)
        213

    All three profiles should have the same product (represent the same lattice volume)

        >>> gsa_profile = GSA(d, n, q, beta, xi, tau)
        >>> cn11_profile = CN11(d, n, q, beta, xi, tau)
        >>> sum([log(x) for x in cn11_profile]
        1296.1852276471009
        >>> sum([log(x) for x in zgsa_profile])
        1296.18522764710
        >>> sum([log(x) for x in gsa_profile])
        1296.18522764710

    Changing xi will change the volume of the lattice

        >>> xi = 2
        >>> gsa_profile = GSA(d, n, q, beta, xi, tau)
        >>> zgsa_profile = ZGSA(d, n, q, beta, xi, tau)
        >>> cn11_profile = CN11(d, n, q, beta, xi, tau)
        >>> sum([log(x) for x in gsa_profile])
        1473.63090587044
        >>> sum([log(x) for x in zgsa_profile])
        1473.63090587044
        >>> sum([log(x) for x in cn11_profile])
        1473.630905870442
    """

    @cached_function
    def ball_log_vol(n):
        return RR((n/2.) * log(pi) - lgamma(n/2. + 1))

    def log_gh(d, logvol=0):
        if d < 49:
            return RR(gh_constant[d] + logvol/d)

        return RR(1./d * (logvol - ball_log_vol(d)))

    def delta(k):
        assert k >= 60
        delta = exp(log_gh(k)/(k-1))
        return RR(delta)

    @cached_function
    def slope(beta):
        if beta<=60:
            return small_slope_t8[beta]
        if beta<=70:
            # interpolate between experimental and asymptotics
            ratio = (70-beta)/10.
            return ratio*small_slope_t8[60]+(1.-ratio)*2*log(delta(70))
        else:
            return 2 * log(delta(beta))

    if not tau:
        L_log = (d - n)*[RR(log(q))] + n * [RR(log(xi))]
        num_q_vec = (d - n)

    else:
        L_log = (d - n - 1)*[RR(log(q))] + n * [RR(log(xi))] + [RR(log(tau))]
        num_q_vec = (d - n - 1)

    slope_ = slope(beta)
    diff = slope(beta)/2.

    for i in range(num_q_vec):
        if diff > (RR(log(q)) - RR(log(xi)))/2.:
            break

        low = (num_q_vec)-i-1
        high = (num_q_vec) + i
        if low >= 0:
            L_log[low] = (RR(log(q)) + RR(log(xi)))/2. + diff

        if high < len(L_log):
            L_log[high] = (RR(log(q)) + RR(log(xi)))/2. - diff

        diff += slope_

    # Output basis profile as squared lengths, not ln(length)
    L = [exp(2 * l_) for l_ in sorted(L_log, reverse=True)]
    return L


def LGSA(d, n, q, beta, xi=1, tau=1, dual=False):
    """
    Reduced lattice shape following the Z-shape Geometric Series Assumption with basis
    rerandomization. Results in BKZ 'forgetting' the q-vectors []_

    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n-1`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: Kannan factor τ.
    :param dual: ignored, since LGSA is self-dual: applying the GSA to the dual is equivalent to
           applying it to the primal.
    :returns: squared Gram-Schmidt norms

    """
    from .reduction import delta as deltaf

    if not tau:
        log_vol = RR((d - n)*log(q, 2) + n*log(xi, 2))
        r_log = d*[RR(log(xi, 2))]
        profile_log_vol = d*RR(log(xi, 2))

    else:
        log_vol = RR((d - n - 1)*log(q, 2) + n*log(xi, 2) + log(tau, 2))
        r_log = (d - 1)*[RR(log(xi, 2))] + [RR(log(tau, 2))]
        profile_log_vol = RR((d - 1)*log(xi, 2) + log(tau, 2))

    slope = -2 * RR(log(deltaf(beta), 2))
    log_vec_len = 0
    for i in range(d - 1, -1, -1):
        log_vec_len -= slope
        profile_log_vol += log_vec_len
        r_log[i] += log_vec_len

        if profile_log_vol > log_vol:
            break

    # Rearrange r to have proper profile shape
    num_gsa_vec = d - i
    r_log = sorted(r_log, reverse=True)

    profile_log_vol = sum(r_log)
    diff = profile_log_vol - log_vol

    for i in range(num_gsa_vec):  # Small shift of the GSA sequence to fix volume
        r_log[i] -= diff / num_gsa_vec

    profile_log_vol = sum(r_log)
    assert abs(profile_log_vol/log_vol - 1) < 1e-6  # Sanity check the volume

    r = [2**(2 * r_) for r_ in r_log]
    return r


def normalize(name):
    if str(name).upper() == "CN11":
        return CN11

    if str(name).upper() == "CN11_GSA":
        return partial(CN11, ignore_qary=True)

    if str(name).upper() == "GSA":
        return GSA

    if str(name).upper() == "ZGSA":
        return ZGSA

    if str(name).upper() == "LGSA":
        return LGSA

    return name


def plot_gso(r, *args, **kwds):
    return line([(i, log(r_, 2) / 2.0) for i, r_ in enumerate(r)], *args, **kwds)
