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
    if tau is None:
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

    if tau is None:
        log_vol = RR(log(q, 2) * (d - n) + log(xi, 2) * n)
    else:
        log_vol = RR(log(q, 2) * (d - n - 1) + log(xi, 2) * n + log(tau, 2))

    delta = deltaf(beta)
    r_log = [(d - 1 - 2 * i) * RR(log(delta, 2)) + log_vol / d for i in range(d)]
    r = [2 ** (2 * r_) for r_ in r_log]
    return r


def ZGSA(d, n, q, beta, xi=1, tau=1, dual=False):
    from math import lgamma
    """
    Reduced lattice Z-shape following the Geometric Series Assumption as specified in
    NTRU fatrigue [DucWoe21]
    :param d: Lattice dimension.
    :param n: The number of `q` vectors is `d-n`.
    :param q: Modulus `q`
    :param beta: Block size β.
    :param xi: Scaling factor ξ for identity part.
    :param tau: ignored, as NTRU represents a homogenous system
    :param dual: ignored, since GSA is self-dual: applying the GSA to the dual is equivalent to
           applying it to the primal.
    :returns: Log Gram-Schmidt norms
    """
    gh_constant = {1:0.00000,2:-0.50511,3:-0.46488,4:-0.39100,5:-0.29759,6:-0.24880,7:-0.21970,8:-0.15748,9:-0.14673,10:-0.07541,11:-0.04870,12:-0.01045,13:0.02298,14:0.04212,15:0.07014,16:0.09205,17:0.12004,18:0.14988,19:0.17351,20:0.18659,21:0.20971,22:0.22728,23:0.24951,24:0.26313,25:0.27662,26:0.29430,27:0.31399,28:0.32494,29:0.34796,30:0.36118,31:0.37531,32:0.39056,33:0.39958,34:0.41473,35:0.42560,36:0.44222,37:0.45396,38:0.46275,39:0.47550,40:0.48889,41:0.50009,42:0.51312,43:0.52463,44:0.52903,45:0.53930,46:0.55289,47:0.56343,48:0.57204,49:0.58184,50:0.58852}
    small_slope_t8 = {2:0.04473,3:0.04472,4:0.04402,5:0.04407,6:0.04334,7:0.04326,8:0.04218,9:0.04237,10:0.04144,11:0.04054,12:0.03961,13:0.03862,14:0.03745,15:0.03673,16:0.03585,17:0.03477,18:0.03378,19:0.03298,20:0.03222,21:0.03155,22:0.03088,23:0.03029,24:0.02999,25:0.02954,26:0.02922,27:0.02891,28:0.02878,29:0.02850,30:0.02827,31:0.02801,32:0.02786,33:0.02761,34:0.02768,35:0.02744,36:0.02728,37:0.02713,38:0.02689,39:0.02678,40:0.02671,41:0.02647,42:0.02634,43:0.02614,44:0.02595,45:0.02583,46:0.02559,47:0.02534,48:0.02514,49:0.02506,50:0.02493,51:0.02475,52:0.02454,53:0.02441,54:0.02427,55:0.02407,56:0.02393,57:0.02371,58:0.02366,59:0.02341,60:0.02332}


    @cached_function
    def ball_log_vol(n):
        return RR((n/2.) * log(pi) - lgamma(n/2. + 1))


    def log_gh(d, logvol=0):
        if d < 49:
            return RR(gh_constant[d] + logvol/d)

        return RR(1./d * (logvol - ball_log_vol(d)))


    def delta(k):
        assert(k>=60)
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

    L_log = (d - n)*[logq] + n * [0]
    slope_ = slope(beta)
    diff = slope(beta)/2.

    for i in range(d//2):
        if diff > logq/2.: break
        try:
            L_log[(d//2)-i-1] = logq/2. + diff
            L_log[(d//2)+i  ] = logq/2. - diff

        except IndexError:
            print(f"d:{d}, n:{n}, d-n: {d-n}, i:{i}, {d-n-i-1}, {d-n+i}, len[l]: {len(L_log)}")

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
