# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using dial attacks.

See :ref:`LWE Dual Attacks` for an introduction what is available.

"""

from functools import partial
from dataclasses import replace

from sage.all import oo, ceil, sqrt, log, cached_function
from .reduction import delta as deltaf
from .util import local_minimum
from .cost import Cost
from .lwe_parameters import LWEParameters
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .io import Logging
from .conf import red_cost_model as red_cost_model_default
from .conf import mitm_opt as mitm_opt_default
from .errors import OutOfBoundsError
from .nd import NoiseDistribution
from .lwe_guess import exhaustive_search, mitm, distinguish


class DualHybrid:
    """
    Estimate cost of solving LWE using dual attacks.
    """

    @staticmethod
    @cached_function
    def dual_reduce(
        delta: float,
        params: LWEParameters,
        zeta: int = 0,
        h1: int = 0,
        rho: float = 1.0,
        log_level=None,
    ):
        """
        Produce new LWE sample using a dual vector on first `n-ζ` coordinates of the secret. The
        length of the dual vector is given by `δ` in root Hermite form and using a possible
        scaling factor, i.e. `|v| = ρ ⋅ δ^d * q^((n-ζ)/d)`.

        :param delta: Length of the vector in root Hermite form
        :param params: LWE parameters
        :param zeta: Dimension ζ ≥ 0 of new LWE instance
        :param h1: Number of non-zero components of the secret of the new LWE instance
        :param rho: Factor introduced by obtaining multiple dual vectors
        :returns: new ``LWEParameters`` and ``m``

        .. note :: This function assumes that the instance is normalized.

        """
        if not 0 <= zeta <= params.n:
            raise OutOfBoundsError(
                f"Splitting dimension {zeta} must be between 0 and n={params.n}."
            )

        # Compute new secret distribution

        if params.Xs.is_sparse:
            h = params.Xs.get_hamming_weight(params.n)
            if not 0 <= h1 <= h:
                raise OutOfBoundsError(f"Splitting weight {h1} must be between 0 and h={h}.")
            # assuming the non-zero entries are uniform
            p = h1 / 2
            red_Xs = NoiseDistribution.SparseTernary(params.n - zeta, h / 2 - p)
            slv_Xs = NoiseDistribution.SparseTernary(zeta, p)

            if h1 == h:
                # no reason to do lattice reduction if we assume
                # that the hw on the reduction part is 0
                return replace(params, Xs=slv_Xs, m=oo), 1
        else:
            # distribution is i.i.d. for each coordinate
            red_Xs = replace(params.Xs, n=params.n - zeta)
            slv_Xs = replace(params.Xs, n=zeta)

        c = red_Xs.stddev * params.q / params.Xe.stddev

        # see if we have optimally many samples (as in [INDOCRYPT:EspJouKha20]) available
        m_ = max(1, ceil(sqrt(red_Xs.n * log(c) / log(delta))) - red_Xs.n)
        m_ = min(params.m, m_)

        # Compute new noise as in [INDOCRYPT:EspJouKha20]
        sigma_ = rho * red_Xs.stddev * delta ** (m_ + red_Xs.n) / c ** (m_ / (m_ + red_Xs.n))
        slv_Xe = NoiseDistribution.DiscreteGaussian(params.q * sigma_)

        slv_params = LWEParameters(
            n=zeta,
            q=params.q,
            Xs=slv_Xs,
            Xe=slv_Xe,
        )

        # The m_ we compute there is the optimal number of samples that we pick from the input LWE
        # instance. We then need to return it because it determines the lattice dimension for the
        # reduction.

        return slv_params, m_

    @staticmethod
    @cached_function
    def cost(
        solver,
        params: LWEParameters,
        beta: int,
        zeta: int = 0,
        h1: int = 0,
        success_probability: float = 0.99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        log_level=None,
    ):
        """
        Computes the cost of the dual hybrid attack that dual reduces the LWE instance and then
        uses the given solver to solve the reduced instance.

        :param solver: Algorithm for solving the reduced instance
        :param params: LWE parameters
        :param beta: Block size used to produce short dual vectors for dual reduction
        :param zeta: Dimension ζ ≥ 0 of new LWE instance
        :param h1: Number of non-zero components of the secret of the new LWE instance
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction
        :param use_lll: Use LLL calls to produce more small vectors

        .. note :: This function assumes that the instance is normalized. It runs no optimization,
            it merely reports costs.

        """
        Logging.log("dual", log_level, f"β={beta}, ζ={zeta}, h1={h1}")

        delta = deltaf(beta)

        # only care about the scaling factor and don't know d yet -> use beta as dummy d
        rho, _, _ = red_cost_model.short_vectors(beta=beta, d=beta)

        params_slv, m_ = DualHybrid.dual_reduce(
            delta, params, zeta, h1, rho, log_level=log_level + 1
        )
        Logging.log("dual", log_level + 1, f"red LWE instance: {repr(params_slv)}")

        cost = solver(params_slv, success_probability)
        Logging.log("dual", log_level + 2, f"solve: {repr(cost)}")

        if cost["rop"] == oo or cost["m"] == oo:
            return replace(cost, beta=beta)

        d = m_ + params.n - zeta
        _, cost_red, _ = red_cost_model.short_vectors(beta, d, cost["m"])
        Logging.log("dual", log_level + 2, f"red: {repr(Cost(rop=cost_red))}")

        cost["rop"] += cost_red
        cost["m"] = m_
        cost["beta"] = beta

        if d < params.n - zeta:
            raise RuntimeError(f"{d} < {params.n - zeta}, {params.n}, {zeta}, {m_}")
        cost["d"] = d

        Logging.log("dual", log_level, f"{repr(cost)}")

        rep = 1
        if params.Xs.is_sparse:
            h = params.Xs.get_hamming_weight(params.n)
            probability = prob_drop(params.n, h, zeta, h1)
            rep = prob_amplify(success_probability, probability)
        # don't need more samples to re-run attack, since we may
        # just guess different components of the secret
        return cost.repeat(times=rep, select={"m": False})

    @staticmethod
    def optimize_blocksize(
        solver,
        params: LWEParameters,
        zeta: int = 0,
        h1: int = 0,
        success_probability: float = 0.99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        log_level=5,
        opt_step=2,
    ):
        """
        Optimizes the cost of the dual hybrid attack over the block size β.

        :param solver: Algorithm for solving the reduced instance
        :param params: LWE parameters
        :param zeta: Dimension ζ ≥ 0 of new LWE instance
        :param h1: Number of non-zero components of the secret of the new LWE instance
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction
        :param use_lll: Use LLL calls to produce more small vectors
        :param opt_step: control robustness of optimizer

        .. note :: This function assumes that the instance is normalized. ζ and h1 are fixed.

        """

        f = partial(
            DualHybrid.cost,
            solver=solver,
            params=params,
            zeta=zeta,
            h1=h1,
            success_probability=success_probability,
            red_cost_model=red_cost_model,
            use_lll=use_lll,
            log_level=log_level,
        )
        # don't have a reliable upper bound for beta
        # we choose n - k arbitrarily and adjust later if
        # necessary
        beta_upper = max(params.n - zeta, 40)
        beta = beta_upper
        while beta == beta_upper:
            beta_upper *= 2
            with local_minimum(40, beta_upper, opt_step) as it:
                for beta in it:
                    it.update(f(beta=beta))
                for beta in it.neighborhood:
                    it.update(f(beta=beta))
                cost = it.y
            beta = cost["beta"]

        cost["zeta"] = zeta
        if params.Xs.is_sparse:
            cost["h1"] = h1
        return cost

    def __call__(
        self,
        solver,
        params: LWEParameters,
        success_probability: float = 0.99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        opt_step=2,
        log_level=1,
    ):
        """
        Optimizes the cost of the dual hybrid attack (using the given solver) over
        all attack parameters: block size β, splitting dimension ζ, and
        splitting weight h1 (in case the secret distribution is sparse). Since
        the cost function for the dual hybrid might only be convex in an approximate
        sense, the parameter ``opt_step`` allows to make the optimization procedure more
        robust against local irregularities (higher value) at the cost of a longer
        running time. In a nutshell, if the cost of the dual hybrid seems suspiciosly
        high, try a larger ``opt_step`` (e.g. 4 or 8).

        :param solver: Algorithm for solving the reduced instance
        :param params: LWE parameters
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction
        :param use_lll: use LLL calls to produce more small vectors [EC:Albrecht17]_
        :param opt_step: control robustness of optimizer

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: Total amount of memory used by solver (in elements mod q).
        - ``red``: Number of word operations in lattice reduction.
        - ``δ``: Root-Hermite factor targeted by lattice reduction.
        - ``β``: BKZ block size.
        - ``ζ``: Number of guessed coordinates.
        - ``h1``: Number of non-zero components among guessed coordinates (if secret distribution is sparse)
        - ``prob``: Probability of success in guessing.
        - ``repetitions``: How often we are required to repeat the attack.
        - ``d``: Lattice dimension.

        - When ζ = 1 this function essentially estimates the dual attack.
        - When ζ > 1 and ``solver`` is ``exhaustive_search`` this function estimates
            the hybrid attack as given in [INDOCRYPT:EspJouKha20]_
        - When ζ > 1 and ``solver`` is ``mitm`` this function estimates the dual MITM
            hybrid attack roughly following [EPRINT:CHHS19]_

        EXAMPLES::

            >>> from estimator import *
            >>> params = LWE.Parameters(n=1024, q = 2**32, Xs=ND.Uniform(0,1), Xe=ND.DiscreteGaussian(3.0))
            >>> LWE.dual(params)
            rop: ≈2^103.4, mem: ≈2^47.9, m: 982, β: 269, d: 2006, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^99.3, mem: ≈2^93.9, m: 944, β: 253, d: 1921, ↻: 1, ζ: 47, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^125.1, mem: ≈2^121.9, m: 1153, k: 115, ↻: 1, β: 351, d: 2039, ζ: 138, tag: dual_mitm_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization="numerical")
            rop: ≈2^125.1, m: 1154, k: 111, mem: ≈2^118.8, ↻: 1, β: 351, d: 2048, ζ: 130, tag: dual_mitm_hybrid

            >>> params = params.updated(Xs=ND.SparseTernary(params.n, 32))
            >>> LWE.dual(params)
            rop: ≈2^100.0, mem: ≈2^45.6, m: 916, β: 256, d: 1940, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^89.8, mem: ≈2^77.5, m: 743, β: 181, d: 1518, ↻: 582, ζ: 249, h1: 8, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^96.4, mem: ≈2^79.6, m: 771, k: 277, ↻: ≈2^15.8, β: 198, d: 1349, ζ: 446, h1: 17...

            >>> params = params.updated(Xs=ND.CenteredBinomial(8))
            >>> LWE.dual(params)
            rop: ≈2^110.5, mem: ≈2^53.0, m: 1113, β: 295, d: 2137, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^109.6, mem: ≈2^99.1, m: 1105, β: 292, d: 2118, ↻: 1, ζ: 11, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^148.3, mem: ≈2^142.0, m: 1416, k: 33, ↻: 1, β: 439, d: 2406, ζ: 34, tag: dual_mitm_hybrid

            >>> params = params.updated(Xs=ND.DiscreteGaussian(3.0))
            >>> LWE.dual(params)
            rop: ≈2^112.2, mem: ≈2^56.0, m: 1150, β: 302, d: 2174, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^111.9, mem: ≈2^103.3, m: 1147, β: 301, d: 2163, ↻: 1, ζ: 8, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^153.1, mem: ≈2^138.7, m: 1478, k: 22, ↻: 1, β: 458, d: 2480, ζ: 22, tag: dual_mitm_hybrid
        """

        Cost.register_impermanent(
            rop=True,
            mem=False,
            red=True,
            beta=False,
            delta=False,
            m=True,
            d=False,
            zeta=False,
        )
        Logging.log("dual", log_level, f"costing LWE instance: {repr(params)}")

        params = params.normalize()

        if params.Xs.is_sparse:
            Cost.register_impermanent(h1=False)

            def _optimize_blocksize(
                solver,
                params: LWEParameters,
                zeta: int = 0,
                success_probability: float = 0.99,
                red_cost_model=red_cost_model_default,
                use_lll=True,
                log_level=None,
            ):
                h = params.Xs.get_hamming_weight(params.n)
                h1_min = max(0, h - (params.n - zeta))
                h1_max = min(zeta, h)
                Logging.log("dual", log_level, f"h1 ∈ [{h1_min},{h1_max}] (zeta={zeta})")
                with local_minimum(h1_min, h1_max, log_level=log_level + 1) as it:
                    for h1 in it:
                        cost = self.optimize_blocksize(
                            h1=h1,
                            solver=solver,
                            params=params,
                            zeta=zeta,
                            success_probability=success_probability,
                            red_cost_model=red_cost_model,
                            use_lll=use_lll,
                            log_level=log_level + 2,
                        )
                        it.update(cost)
                    return it.y

        else:
            _optimize_blocksize = self.optimize_blocksize

        f = partial(
            _optimize_blocksize,
            solver=solver,
            params=params,
            success_probability=success_probability,
            red_cost_model=red_cost_model,
            use_lll=use_lll,
            log_level=log_level + 1,
        )

        with local_minimum(1, params.n - 1, opt_step) as it:
            for zeta in it:
                it.update(f(zeta=zeta))
            for zeta in it.neighborhood:
                it.update(f(zeta=zeta))
            cost = it.y

        cost["problem"] = params
        return cost


DH = DualHybrid()


def dual(
    params: LWEParameters,
    success_probability: float = 0.99,
    red_cost_model=red_cost_model_default,
    use_lll=True,
):
    """
    Dual hybrid attack as in [PQCBook:MicReg09]_.

    :param params: LWE parameters.
    :param success_probability: The success probability to target.
    :param red_cost_model: How to cost lattice reduction.
    :param use_lll: use LLL calls to produce more small vectors [EC:Albrecht17]_.

    The returned cost dictionary has the following entries:

    - ``rop``: Total number of word operations (≈ CPU cycles).
    - ``mem``: Total amount of memory used by solver (in elements mod q).
    - ``red``: Number of word operations in lattice reduction.
    - ``δ``: Root-Hermite factor targeted by lattice reduction.
    - ``β``: BKZ block size.
    - ``prob``: Probability of success in guessing.
    - ``repetitions``: How often we are required to repeat the attack.
    - ``d``: Lattice dimension.

    """
    Cost.register_impermanent(
        rop=True,
        mem=False,
        red=True,
        beta=False,
        delta=False,
        m=True,
        d=False,
    )

    ret = DH.optimize_blocksize(
        solver=distinguish,
        params=params,
        zeta=0,
        h1=0,
        success_probability=success_probability,
        red_cost_model=red_cost_model,
        use_lll=use_lll,
        log_level=1,
    )
    del ret["zeta"]
    if hasattr(ret, "h1"):
        del ret["h1"]
    ret["tag"] = "dual"
    return ret


def dual_hybrid(
    params: LWEParameters,
    success_probability: float = 0.99,
    red_cost_model=red_cost_model_default,
    use_lll=True,
    mitm_optimization=False,
    opt_step=2,
):
    """
    Dual hybrid attack from [INDOCRYPT:EspJouKha20]_.

    :param params: LWE parameters.
    :param success_probability: The success probability to target.
    :param red_cost_model: How to cost lattice reduction.
    :param use_lll: Use LLL calls to produce more small vectors [EC:Albrecht17]_.
    :param mitm_optimization: One of "analytical" or "numerical". If ``True`` a default from the
           ``conf`` module is picked, ``False`` disables MITM.
    :param opt_step: Control robustness of optimizer.

    The returned cost dictionary has the following entries:

    - ``rop``: Total number of word operations (≈ CPU cycles).
    - ``mem``: Total amount of memory used by solver (in elements mod q).
    - ``red``: Number of word operations in lattice reduction.
    - ``δ``: Root-Hermite factor targeted by lattice reduction.
    - ``β``: BKZ block size.
    - ``ζ``: Number of guessed coordinates.
    - ``h1``: Number of non-zero components among guessed coordinates (if secret distribution is sparse)
    - ``prob``: Probability of success in guessing.
    - ``repetitions``: How often we are required to repeat the attack.
    - ``d``: Lattice dimension.
    """

    if mitm_optimization is True:
        mitm_optimization = mitm_opt_default

    if mitm_optimization:
        solver = partial(mitm, optimization=mitm_optimization)
    else:
        solver = exhaustive_search

    ret = DH(
        solver=solver,
        params=params,
        success_probability=success_probability,
        red_cost_model=red_cost_model,
        use_lll=use_lll,
        opt_step=opt_step,
    )
    if mitm_optimization:
        ret["tag"] = "dual_mitm_hybrid"
    else:
        ret["tag"] = "dual_hybrid"
    return ret
