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


def qary_simulator(f, d, n, q, beta, xi=1, tau=1, dual=False):
    """
    Reduced lattice shape calling ``f``.

    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n-1`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: Kannan factor τ.
    :param dual: perform reduction on the dual.

    """
    if not tau:
        r = [q**2] * (d - n) + [xi**2] * n
    else:
        r = [q**2] * (d - n - 1) + [xi**2] * n + [tau**2]

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


def CN11(d, n, q, beta, xi=1, tau=1, dual=False):
    """
    Reduced lattice shape using simulator from [AC:CheNgu11]_

    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n-1`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: Kannan factor τ.
    :param dual: perform reduction on the dual.
    :returns: squared Gram-Schmidt norms

    """

    from fpylll import BKZ
    from fpylll.tools.bkz_simulator import simulate

    def f(r, beta):
        return simulate(r, BKZ.EasyParam(beta))[0]

    return qary_simulator(f=f, d=d, n=n, q=q, beta=beta, xi=xi, tau=tau, dual=dual)


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
    NTRU fatrigue [DucWoe21]
    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param dual: ignored, since GSA is self-dual: applying the GSA to the dual is equivalent to
           applying it to the primal.
    :returns: Squared Gram-Schmidt norms
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

    # Scale down q by xi, instead of scaling Identity part up by xi.
    logq = RR(log(q) - log(xi))

    if not tau:
        L_log = (d - n)*[logq] + n * [0]
    else:
        L_log = (d - n)*[logq] + n * [0] + [RR(log(tau))]

    slope_ = slope(beta)
    diff = slope(beta)/2.

    for i in range(d-n):
        if diff > logq/2.:
            break

        low = (d - n)-i-1
        high = (d - n) + i
        if low >= 0:
            L_log[low] = logq/2. + diff

        if high < len(L_log):
            L_log[high] = logq/2. - diff

        diff += slope_

    # Scale basis profile_using xi, to induce scaling on Identity part
    scale_factor = (logq + RR(log(xi)))/logq
    L_log = [scale_factor*l_ for l_ in L_log]

    # Output basis profile as squared lengths, not ln(length)
    L = [exp(2 * l_) for l_ in L_log]
    return L


def normalize(name):
    if str(name).upper() == "CN11":
        return CN11

    if str(name).upper() == "GSA":
        return GSA

    if str(name).upper() == "ZGSA":
        return ZGSA

    return name


def plot_gso(r, *args, **kwds):
    return line([(i, log(r_, 2) / 2.0) for i, r_ in enumerate(r)], *args, **kwds)
