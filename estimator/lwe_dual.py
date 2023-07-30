# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using dial attacks.

See :ref:`LWE Dual Attacks` for an introduction what is available.

"""

from functools import partial
from dataclasses import replace

from sage.rings.all import QQ
from sage.functions.log import exp
from sage.functions.log import log
from sage.functions.other import ceil
from sage.misc.cachefunc import cached_function
from sage.misc.functional import sqrt
from sage.rings.infinity import infinity as oo
from sage.rings.real_mpfr import RR
from sage.symbolic.constants import pi

from .reduction import delta as deltaf
from .util import local_minimum
from .cost import Cost
from .lwe_parameters import LWEParameters
from .prob import drop as prob_drop, amplify as prob_amplify
from .io import Logging
from .conf import (red_cost_model as red_cost_model_default,
                   mitm_opt as mitm_opt_default)
from .errors import OutOfBoundsError, InsufficientSamplesError
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
        t: int = 0,
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

        # apply the [AC:GuoJoh21] technique, m_ not optimal anymore?
        d = m_ + red_Xs.n
        rho /= 2 ** (t / d)

        # Compute new noise as in [INDOCRYPT:EspJouKha20]
        # ~ sigma_ = rho * red_Xs.stddev * delta ** (m_ + red_Xs.n) / c ** (m_ / (m_ + red_Xs.n))
        sigma_ = rho * red_Xs.stddev * delta**d / c ** (m_ / d)
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
        t: int = 0,
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

        # only care about the scaling factor and don't know d yet -> use 2 * beta as dummy d
        rho = red_cost_model.short_vectors(beta=beta, d=2 * beta)[0]

        params_slv, m_ = DualHybrid.dual_reduce(
            delta, params, zeta, h1, rho, t, log_level=log_level + 1
        )
        Logging.log("dual", log_level + 1, f"red LWE instance: {repr(params_slv)}")

        if t:
            cost = DualHybrid.fft_solver(params_slv, success_probability, t)
        else:
            cost = solver(params_slv, success_probability)

        Logging.log("dual", log_level + 2, f"solve: {cost!r}")

        if cost["rop"] == oo or cost["m"] == oo:
            cost["beta"] = beta
            return cost

        d = m_ + params.n - zeta
        cost_red = red_cost_model.short_vectors(beta, d, cost["m"])[1]
        Logging.log("dual", log_level + 2, f"red: {Cost(rop=cost_red)!r}")

        cost["rop"] += cost_red
        cost["m"] = m_
        cost["beta"] = beta
        if t:
            cost["t"] = t

        if d < params.n - zeta:
            raise RuntimeError(f"{d} < {params.n - zeta}, {params.n}, {zeta}, {m_}")
        cost["d"] = d

        Logging.log("dual", log_level, f"{repr(cost)}")

        rep = 1
        if params.Xs.is_sparse:
            h = params.Xs.get_hamming_weight(params.n)
            probability = RR(prob_drop(params.n, h, zeta, h1))
            rep = prob_amplify(success_probability, probability)
        # don't need more samples to re-run attack, since we may
        # just guess different components of the secret
        return cost.repeat(times=rep, select={"m": False})

    @staticmethod
    def fft_solver(params, success_probability, t=0):
        """
        Estimate cost of solving LWE via the FFT distinguisher from [AC:GuoJoh21]_.

        :param params: LWE parameters
        :param success_probability: the targeted success probability
        :param t: the number of secret coordinates to guess mod 2.
            For t=0 this is similar to lwe_guess.ExhaustiveSearch.
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: memory requirement in integers mod q.
        - ``m``: Required number of samples to distinguish the correct solution with high probability.

        .. note :: The parameter t only makes sense in the context of the dual attack,
            which is why this function is here and not in the lwe_guess module.
        """

        # there are two stages: enumeration and distinguishing, so we split up the success_probability
        probability = sqrt(success_probability)

        try:
            size = params.Xs.support_size(n=params.n, fraction=probability)
            size_fft = 2**t
        except NotImplementedError:
            # not achieving required probability with search space
            # given our settings that means the search space is huge
            # so we approximate the cost with oo
            return Cost(rop=oo, mem=oo, m=1)

        sigma = params.Xe.stddev / params.q
        m_required = RR(
            4
            * exp(4 * pi * pi * sigma * sigma)
            * (log(size_fft * size) - log(log(1 / probability)))
        )

        if params.m < m_required:
            raise InsufficientSamplesError(
                f"Exhaustive search: Need {m_required} samples but only {params.m} available."
            )
        else:
            m = m_required

        cost = size * (m + t * size_fft)

        return Cost(rop=cost, mem=cost, m=m)

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
        opt_step=8,
        fft=False,
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
        :param fft: use the FFT distinguisher from [AC:GuoJoh21]_

        .. note :: This function assumes that the instance is normalized. ζ and h1 are fixed.

        """

        f_t = partial(
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

        if fft is True:

            def f(beta):
                with local_minimum(0, params.n - zeta) as it:
                    for t in it:
                        it.update(f_t(beta=beta, t=t))
                    return it.y

        else:
            f = f_t

        # don't have a reliable upper bound for beta
        # we choose n - k arbitrarily and adjust later if
        # necessary
        beta_upper = min(max(params.n - zeta, 40), 1024)
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
        opt_step=8,
        log_level=1,
        fft=False,
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
        :param fft: use the FFT distinguisher from [AC:GuoJoh21]_. (ignored for sparse secrets)

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: Total amount of memory used by solver (in elements mod q).
        - ``red``: Number of word operations in lattice reduction.
        - ``β``: BKZ block size.
        - ``ζ``: Number of guessed coordinates.
        - ``h1``: Number of non-zero components among guessed coordinates (if secret distribution is sparse)
        - ``prob``: Probability of success in guessing.
        - ``repetitions``: How often we are required to repeat the attack.
        - ``d``: Lattice dimension.
        - ``t``: Number of secrets to guess mod 2 (only if ``fft`` is ``True``)

        - When ζ = 1 this function essentially estimates the dual attack.
        - When ζ > 1 and ``solver`` is ``exhaustive_search`` this function estimates
            the hybrid attack as given in [INDOCRYPT:EspJouKha20]_
        - When ζ > 1 and ``solver`` is ``mitm`` this function estimates the dual MITM
            hybrid attack roughly following [EPRINT:CHHS19]_

        EXAMPLES::

            >>> from estimator import *
            >>> params = LWE.Parameters(n=1024, q = 2**32, Xs=ND.Uniform(0,1), Xe=ND.DiscreteGaussian(3.0))
            >>> LWE.dual(params)
            rop: ≈2^107.0, mem: ≈2^58.0, m: 970, β: 264, d: 1994, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^103.2, mem: ≈2^97.4, m: 937, β: 250, d: 1919, ↻: 1, ζ: 42, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^130.1, mem: ≈2^127.0, m: 1144, k: 120, ↻: 1, β: 347, d: 2024, ζ: 144, tag: dual_mitm_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization="numerical")
            rop: ≈2^129.0, m: 1145, k: 1, mem: ≈2^131.0, ↻: 1, β: 346, d: 2044, ζ: 125, tag: dual_mitm_hybrid

            >>> params = params.updated(Xs=ND.SparseTernary(params.n, 32))
            >>> LWE.dual(params)
            rop: ≈2^103.4, mem: ≈2^55.4, m: 904, β: 251, d: 1928, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^92.1, mem: ≈2^78.2, m: 716, β: 170, d: 1464, ↻: 1989, ζ: 276, h1: 8, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^98.2, mem: ≈2^78.6, m: 728, k: 292, ↻: ≈2^18.7, β: 180, d: 1267, ζ: 485, h1: 17, tag: ...

            >>> params = params.updated(Xs=ND.CenteredBinomial(8))
            >>> LWE.dual(params)
            rop: ≈2^114.5, mem: ≈2^61.0, m: 1103, β: 291, d: 2127, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^113.6, mem: ≈2^103.5, m: 1096, β: 288, d: 2110, ↻: 1, ζ: 10, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^155.5, mem: ≈2^146.2, m: 1414, k: 34, ↻: 1, β: 438, d: 2404, ζ: 34, tag: dual_mitm_hybrid

            >>> params = params.updated(Xs=ND.DiscreteGaussian(3.0))
            >>> LWE.dual(params)
            rop: ≈2^116.5, mem: ≈2^64.0, m: 1140, β: 298, d: 2164, ↻: 1, tag: dual
            >>> LWE.dual_hybrid(params)
            rop: ≈2^116.2, mem: ≈2^100.4, m: 1137, β: 297, d: 2155, ↻: 1, ζ: 6, tag: dual_hybrid
            >>> LWE.dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^160.7, mem: ≈2^156.8, m: 1473, k: 25, ↻: 1, β: 456, d: 2472, ζ: 25, tag: dual_mitm_hybrid

            >>> LWE.dual_hybrid(schemes.NTRUHPS2048509Enc)
            rop: ≈2^131.7, mem: ≈2^128.5, m: 436, β: 358, d: 906, ↻: 1, ζ: 38, tag: dual_hybrid

            >>> LWE.dual(schemes.CHHS_4096_67)
            rop: ≈2^206.9, mem: ≈2^126.0, m: ≈2^11.8, β: 616, d: 7779, ↻: 1, tag: dual

            >>> LWE.dual_hybrid(schemes.Kyber512, red_cost_model=RC.GJ21, fft=True)
            rop: ≈2^149.6, mem: ≈2^145.7, m: 510, β: 399, t: 76, d: 1000, ↻: 1, ζ: 22, tag: dual_hybrid
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
            t=False,
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
                fft=False,
            ):
                h = params.Xs.get_hamming_weight(params.n)
                h1_min = max(0, h - (params.n - zeta))
                h1_max = min(zeta, h)
                if h1_min == h1_max:
                    h1_max = h1_min + 1
                Logging.log("dual", log_level, f"h1 ∈ [{h1_min},{h1_max}] (zeta={zeta})")
                with local_minimum(h1_min, h1_max, log_level=log_level + 1) as it:
                    for h1 in it:
                        # ignoring fft on purpose for sparse secrets
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
            fft=fft,
        )

        with local_minimum(1, params.n - 1, opt_step) as it:
            for zeta in it:
                it.update(f(zeta=zeta))
            for zeta in it.neighborhood:
                it.update(f(zeta=zeta))
            cost = it.y

        cost["problem"] = params
        return cost.sanity_check()


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
    if "h1" in ret:
        del ret["h1"]
    ret["tag"] = "dual"
    return ret


def dual_hybrid(
    params: LWEParameters,
    success_probability: float = 0.99,
    red_cost_model=red_cost_model_default,
    use_lll=True,
    mitm_optimization=False,
    opt_step=8,
    fft=False,
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
    :param fft: use the FFT distinguisher from [AC:GuoJoh21]_. (ignored for sparse secrets)

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
    - ``t``: Number of secrets to guess mod 2 (only if ``fft`` is ``True``)
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
        fft=fft,
    )
    if mitm_optimization:
        ret["tag"] = "dual_mitm_hybrid"
    else:
        ret["tag"] = "dual_hybrid"
    return ret
