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

from sage.functions.log import log
from sage.plot.line import line
from sage.rings.real_mpfr import RR

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


def normalize(name):
    if str(name).upper() == "CN11":
        return CN11

    if str(name).upper() == "GSA":
        return GSA

    return name


def plot_gso(r, *args, **kwds):
    return line([(i, log(r_, 2) / 2.0) for i, r_ in enumerate(r)], *args, **kwds)
