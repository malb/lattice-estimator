# -*- coding: utf-8 -*-
"""
Default values.
"""

from .reduction import RC
from .simulator import GSA
from sage.all import exp

red_cost_model = RC.MATZOV
red_cost_model_classical_poly_space = RC.ABLR21
red_shape_model = "gsa"
red_simulator = GSA
mitm_opt = "analytical"
ntru_fatigue_lb = lambda n: int((n**2.484)/exp(6))
ntru_fatigue_ub = lambda n: int((n**2.484)/exp(5))
