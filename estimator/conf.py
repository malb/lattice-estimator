# -*- coding: utf-8 -*-
"""
Default values.
"""

from .reduction import Kyber, ABLR21
from .simulator import GSA

red_cost_model_default = Kyber
red_cost_model_classical_poly_space = ABLR21
red_shape_model_default = "gsa"
red_simulator_default = GSA

default_mitm_opt = "analytical"
