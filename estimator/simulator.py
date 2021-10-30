# -*- coding: utf-8 -*-
from sage.all import RR, log


class Simulator:
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

    @classmethod
    def CN11(cls, d, n, q, beta, xi=1, tau=1):
        from fpylll import BKZ
        from fpylll.tools.bkz_simulator import simulate

        if tau is not None:
            r = [q ** 2] * (d - n - 1) + [xi ** 2] * n + [tau ** 2]
        else:
            r = [q ** 2] * (d - n) + [xi ** 2] * n

        return simulate(r, BKZ.EasyParam(beta))[0]

    @classmethod
    def GSA(cls, d, n, q, beta, xi=1, tau=1):
        from .reduction import BKZ

        if tau is not None:
            log_vol = RR(log(q, 2) * (d - n - 1) + log(xi, 2) * n + log(tau, 2))
        else:
            log_vol = RR(log(q, 2) * (d - n) + log(xi, 2) * n)

        delta = BKZ.delta(beta)
        r_log = [(d - 1 - 2 * i) * RR(log(delta, 2)) + log_vol / d for i in range(d)]
        r = [2 ** (2 * r_) for r_ in r_log]
        return r


def plot_gso(r, *args, **kwds):
    from sage.all import line

    return line([(i, log(r_, 2) / 2.0) for i, r_ in enumerate(r)], *args, **kwds)
