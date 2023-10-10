# -*- coding: utf-8 -*-
"""
Default values.
"""

from .simulator import GSA
from .reduction import RC
from sage.all import exp

red_cost_model = RC.MATZOV
red_cost_model_classical_poly_space = RC.ABLR21
red_shape_model = "gsa"
red_simulator = GSA
mitm_opt = "analytical"
max_n_cache = 10000


def ntru_fatigue_lb(n):
    return int((n**2.484)/exp(6))
