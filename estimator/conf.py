# -*- coding: utf-8 -*-
"""
Default values.
"""

from .reduction import RC
from .simulator import GSA

red_cost_model = RC.Kyber
red_cost_model_classical_poly_space = RC.ABLR21
red_shape_model = "gsa"
red_simulator = GSA
mitm_opt = "analytical"
