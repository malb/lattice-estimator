# -*- coding: utf-8 -*-
"""
Default values.
"""

from .simulator import GSA
from .reduction import RC
from sage.all import exp

"""
Default models used to evaluate the cost and shape of lattice reduction.
This influences the concrete estimated cost of attacks.
"""
red_cost_model = RC.MATZOV
red_cost_model_classical_poly_space = RC.ABLR21
red_shape_model = "gsa"
red_simulator = GSA

# Upper bound on blocksize, for security levels of at most 512 bits.
# This value is selected as RC.ADPS16(1754, 1754) ~ 2^(512)
max_beta = 1754

mitm_opt = "analytical"
max_n_cache = 10000


def ntru_fatigue_lb(n):
    return int((n**2.484)/exp(6))
