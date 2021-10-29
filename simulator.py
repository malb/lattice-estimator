# -*- coding: utf-8 -*-


def CN11Simulator(d, n, q, beta, scale=1, kannan_coeff=1):
    r = [q ** 2] * (d - n - 1) + [scale ** 2] * n + [kannan_coeff ** 2]
    from fpylll import BKZ
    from fpylll.tools.bkz_simulator import simulate

    return simulate(r, BKZ.EasyParam(beta))[0]
