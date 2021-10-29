# -*- coding: utf-8 -*-
from sage.all import RR, log


class Simulator:
    @classmethod
    def CN11(cls, d, n, q, beta, scale=1, kannan_coeff=1):
        from fpylll import BKZ
        from fpylll.tools.bkz_simulator import simulate

        r = [q ** 2] * (d - n - 1) + [scale ** 2] * n + [kannan_coeff ** 2]

        return simulate(r, BKZ.EasyParam(beta))[0]

    @classmethod
    def GSA(cls, d, n, q, beta, scale=1, kannan_coeff=1):
        from .reduction import BKZ

        log_vol = RR(log(q, 2) * (d - n - 1) + log(scale, 2) * n + log(kannan_coeff, 2))
        delta = BKZ.delta(beta)
        r_log = [(d - 1 - 2 * i) * RR(log(delta, 2)) + log_vol / d for i in range(d)]
        r = [2 ** (2 * r_) for r_ in r_log]
        return r
