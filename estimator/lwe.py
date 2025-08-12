# -*- coding: utf-8 -*-
"""
High-level LWE interface
"""

from functools import partial
from sage.all import oo

from .lwe_primal import primal_usvp, primal_bdd, primal_hybrid
from .lwe_bkw import coded_bkw
from .lwe_guess import exhaustive_search, mitm, distinguish, guess_composition  # noqa
from .lwe_dual import dual
from .lwe_dual import matzov as dual_hybrid
from .gb import arora_gb  # noqa
from .lwe_parameters import LWEParameters as Parameters  # noqa
from .conf import (
    red_cost_model as red_cost_model_default,
    red_shape_model as red_shape_model_default,
)
from .util import batch_estimate, f_name
from .reduction import RC
from .io import Logging


class Estimate:

    def rough(self, params, jobs=1, catch_exceptions=True, quiet=False):
        """
        This function makes the following (non-default) somewhat routine assumptions to evaluate the cost of lattice
        reduction, and to provide comparable numbers with most of the literature:

        - The GSA holds.
        - The Core-SVP model holds.

        Provided numbers are notably not directly comparable with the rest of our API, when using the default cost
        models.

        This function furthermore assumes the following heuristics:

        - The primal hybrid attack only applies to sparse secrets.
        - The dual hybrid MITM attack only applies to sparse secrets.
        - Arora-GB only applies to bounded noise with at least `n^2` samples.
        - BKW is not competitive.

        :param params: LWE parameters.
        :param jobs: Use multiple threads in parallel.
        :param catch_exceptions: When an estimate fails, just print a warning.
        :param quiet: suppress printing

        EXAMPLE ::

            >>> from estimator import *
            >>> _ = LWE.estimate.rough(schemes.Kyber512)
            usvp                 :: rop: ≈2^118.6, red: ≈2^118.6, δ: 1.003941, β: 406, d: 998, tag: usvp
            dual_hybrid          :: rop: ≈2^115.5, red: ≈2^115.3, guess: ≈2^112.3, β: 395, p: 5, ζ: 0, t: 40, β': 395...

            >>> _ = LWE.estimate.rough(schemes.Kyber512, quiet=True)

        """
        params = params.normalize()

        algorithms = {}

        algorithms["usvp"] = partial(primal_usvp, red_cost_model=RC.ADPS16, red_shape_model="gsa")
        algorithms["dual_hybrid"] = partial(dual_hybrid, red_cost_model=RC.ADPS16)

        if params.m > params.n**2 and params.Xe.is_bounded:
            if params.Xs.is_sparse:
                algorithms["arora-gb"] = guess_composition(arora_gb.cost_bounded)
            else:
                algorithms["arora-gb"] = arora_gb.cost_bounded

        res_raw = batch_estimate(
            params, algorithms.values(), log_level=1, jobs=jobs, catch_exceptions=catch_exceptions
        )
        res_raw = res_raw[params]
        res = {
            algorithm: v
            for algorithm, attack in algorithms.items()
            for k, v in res_raw.items()
            if f_name(attack) == k
        }

        for algorithm in algorithms:
            if algorithm not in res:
                continue
            result = res[algorithm]
            if result["rop"] != oo:
                Logging.print("estimator", int(quiet), f"{algorithm:20s} :: {result!r}")

        return res

    def __call__(
        self,
        params,
        red_cost_model=red_cost_model_default,
        red_shape_model=red_shape_model_default,
        deny_list=tuple(),
        add_list=tuple(),
        jobs=1,
        catch_exceptions=True,
        quiet=False,
    ):
        """
        Run all estimates, based on the default cost and shape models for lattice reduction.

        :param params: LWE parameters.
        :param red_cost_model: How to cost lattice reduction.
        :param red_shape_model: How to model the shape of a reduced basis (applies to primal attacks)
        :param deny_list: skip these algorithms
        :param add_list: add these ``(name, function)`` pairs to the list of algorithms to estimate.a
        :param jobs: Use multiple threads in parallel.
        :param catch_exceptions: When an estimate fails, just print a warning.
        :param quiet: suppress printing

        EXAMPLE ::

            >>> from estimator import *
            >>> _ = LWE.estimate(schemes.Kyber512)
            bkw                  :: rop: ≈2^178.8, m: ≈2^166.8, mem: ≈2^167.8, b: 14, t1: 0, t2: 16, ℓ: 13, #cod: 448...
            usvp                 :: rop: ≈2^143.8, red: ≈2^143.8, δ: 1.003941, β: 406, d: 998, tag: usvp
            bdd                  :: rop: ≈2^140.2, red: ≈2^139.4, svp: ≈2^139.0, β: 390, η: 421, d: 1007, tag: bdd
            dual                 :: rop: ≈2^149.9, mem: ≈2^97.1, m: 512, β: 424, d: 1024, ↻: 1, tag: dual
            dual_hybrid          :: rop: ≈2^139.7, red: ≈2^139.5, guess: ≈2^135.9, β: 387, p: 5, ζ: 0, t: 50, β': 391...

            >>> _ = LWE.estimate(schemes.Kyber512, quiet=True)

        """
        params = params.normalize()

        algorithms = {}

        algorithms["arora-gb"] = guess_composition(arora_gb)
        algorithms["bkw"] = coded_bkw

        algorithms["usvp"] = partial(
            primal_usvp, red_cost_model=red_cost_model, red_shape_model=red_shape_model
        )
        algorithms["bdd"] = partial(
            primal_bdd, red_cost_model=red_cost_model, red_shape_model=red_shape_model
        )
        algorithms["bdd_hybrid"] = partial(
            primal_hybrid,
            mitm=False,
            babai=False,
            red_cost_model=red_cost_model,
            red_shape_model=red_shape_model,
        )
        # we ignore the case of mitm=True babai=False for now, due to it being overly-optimistic
        algorithms["bdd_mitm_hybrid"] = partial(
            primal_hybrid,
            mitm=True,
            babai=True,
            red_cost_model=red_cost_model,
            red_shape_model=red_shape_model,
        )
        algorithms["dual"] = partial(dual, red_cost_model=red_cost_model)
        algorithms["dual_hybrid"] = partial(dual_hybrid, red_cost_model=red_cost_model)

        algorithms = {k: v for k, v in algorithms.items() if k not in deny_list}
        algorithms.update(add_list)

        res_raw = batch_estimate(
            params, algorithms.values(), log_level=1, jobs=jobs, catch_exceptions=catch_exceptions
        )
        res_raw = res_raw[params]
        res = {
            algorithm: v
            for algorithm, attack in algorithms.items()
            for k, v in res_raw.items()
            if f_name(attack) == k
        }

        for algorithm in algorithms:
            if algorithm not in res:
                continue
            result = res[algorithm]
            if result["rop"] == oo:
                continue
            if algorithm == "bdd_hybrid" and res["bdd"]["rop"] <= result["rop"]:
                continue
            if algorithm == "bdd_mitm_hybrid" and res["bdd_hybrid"]["rop"] <= result["rop"]:
                continue
            if algorithm == "dual_mitm_hybrid" and res["dual_hybrid"]["rop"] < result["rop"]:
                continue
            Logging.print("estimator", int(quiet), f"{algorithm:20s} :: {result!r}")

        return res


estimate = Estimate()
