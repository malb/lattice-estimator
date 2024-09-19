# -*- coding: utf-8 -*-
"""
Estimate cost of solving LWE using dial attacks.

See :ref:`LWE Dual Attacks` for an introduction what is available.

"""

from functools import partial
from dataclasses import replace

from sage.all import oo, ceil, sqrt, log, cached_function, RR, exp, pi, e, coth, tanh

from .reduction import delta as deltaf
from .util import local_minimum, early_abort_range
from .cost import Cost
from .lwe_parameters import LWEParameters
from .prob import drop as prob_drop, amplify as prob_amplify
from .io import Logging
from .conf import red_cost_model as red_cost_model_default, mitm_opt as mitm_opt_default
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
            h = params.Xs.hamming_weight
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
        cost["beta"] = beta

        if cost["rop"] == oo or cost["m"] == oo:
            return cost

        d = m_ + params.n - zeta
        _, cost_red, N, sieve_dim = red_cost_model.short_vectors(beta, d, cost["m"])
        Logging.log("dual", log_level + 2, f"red: {Cost(rop=cost_red)!r}")

        # Add the runtime cost of sieving in dimension `sieve_dim` possibly multiple times.
        cost["rop"] += cost_red

        # Add the memory cost of storing the `N` dual vectors, using `sieve_dim` many coefficients
        # (mod q) to represent them. Note that short dual vectors may actually be described by less
        # bits because its coefficients are generally small, so this is really an upper bound here.
        cost["mem"] += sieve_dim * N
        cost["m"] = m_

        if d < params.n - zeta:
            raise RuntimeError(f"{d} < {params.n - zeta}, {params.n}, {zeta}, {m_}")
        cost["d"] = d

        Logging.log("dual", log_level, f"{repr(cost)}")

        rep = 1
        if params.Xs.is_sparse:
            h = params.Xs.hamming_weight
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
        - ``t``: the number of secret coordinates to guess mod 2.

        .. note :: The parameter t only makes sense in the context of the dual attack,
            which is why this function is here and not in the lwe_guess module.
        """

        # there are two stages: enumeration and distinguishing, so we split up the success_probability
        probability = sqrt(success_probability)

        try:
            size = params.Xs.support_size(probability)
            size_fft = 2**t
        except NotImplementedError:
            # not achieving required probability with search space
            # given our settings that means the search space is huge
            # so we approximate the cost with oo
            return Cost(rop=oo, mem=oo, m=1)

        sigma = params.Xe.stddev / params.q

        # Here, assume the Independence Heuristic, cf. [ia.cr/2023/302].
        # The minimal number of short dual vectors that is required to distinguish the correct
        # guess with probability at least `probability`:
        m_required = RR(
            4
            * exp(4 * pi * pi * sigma * sigma)
            * (log(size_fft * size) - log(log(1 / probability)))
        )

        if params.m < m_required:
            raise InsufficientSamplesError(
                f"Exhaustive search: Need {m_required} samples but only {params.m} available."
            )

        # Running a fast Walsh--Hadamard transform takes time proportional to t 2^t.
        runtime_cost = size * (t * size_fft)
        # Add the cost of updating the FFT tables for all of the enumeration targets.
        # Use "Efficient Updating of the FFT Input", [MATZOV, §5.4]:
        runtime_cost += size * (4 * m_required)

        # This is the number of entries the table should have. Note that it should support
        # (floating point) numbers in the range [-N, N], if ``N`` is the number of dual vectors.
        # However 32-bit floats are good enough in practice.
        memory_cost = size_fft

        return Cost(rop=runtime_cost, mem=memory_cost, m=m_required, t=t)

    @staticmethod
    def optimize_blocksize(
        solver,
        params: LWEParameters,
        zeta: int = 0,
        h1: int = 0,
        success_probability: float = 0.99,
        red_cost_model=red_cost_model_default,
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
        running time. In a nutshell, if the cost of the dual hybrid seems suspiciously
        high, try a larger ``opt_step`` (e.g. 4 or 8).

        :param solver: Algorithm for solving the reduced instance
        :param params: LWE parameters
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction
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
            >>> from estimator.lwe_dual import dual_hybrid
            >>> params = LWE.Parameters(n=1024, q = 2**32, Xs=ND.Uniform(0,1), Xe=ND.DiscreteGaussian(3.0))
            >>> LWE.dual(params)
            rop: ≈2^107.0, mem: ≈2^66.4, m: 970, β: 264, d: 1994, ↻: 1, tag: dual
            >>> dual_hybrid(params)
            rop: ≈2^103.2, mem: ≈2^97.4, m: 937, β: 250, d: 1919, ↻: 1, ζ: 42, tag: dual_hybrid
            >>> dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^130.1, mem: ≈2^127.0, m: 1144, k: 120, ↻: 1, β: 347, d: 2024, ζ: 144, tag: dual_mitm_hybrid
            >>> dual_hybrid(params, mitm_optimization="numerical")
            rop: ≈2^129.0, m: 1145, k: 1, mem: ≈2^131.0, ↻: 1, β: 346, d: 2044, ζ: 125, tag: dual_mitm_hybrid

            >>> params = params.updated(Xs=ND.SparseTernary(params.n, 32))
            >>> LWE.dual(params)
            rop: ≈2^103.4, mem: ≈2^63.9, m: 904, β: 251, d: 1928, ↻: 1, tag: dual
            >>> dual_hybrid(params)
            rop: ≈2^92.1, mem: ≈2^78.2, m: 716, β: 170, d: 1464, ↻: 1989, ζ: 276, h1: 8, tag: dual_hybrid
            >>> dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^98.2, mem: ≈2^78.6, m: 728, k: 292, ↻: ≈2^18.7, β: 180, d: 1267, ζ: 485, h1: 17, tag: ...

            >>> params = params.updated(Xs=ND.CenteredBinomial(8))
            >>> LWE.dual(params)
            rop: ≈2^114.5, mem: ≈2^71.8, m: 1103, β: 291, d: 2127, ↻: 1, tag: dual
            >>> dual_hybrid(params)
            rop: ≈2^113.6, mem: ≈2^103.5, m: 1096, β: 288, d: 2110, ↻: 1, ζ: 10, tag: dual_hybrid
            >>> dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^155.5, mem: ≈2^146.2, m: 1414, k: 34, ↻: 1, β: 438, d: 2404, ζ: 34, tag: dual_mitm_hybrid

            >>> params = params.updated(Xs=ND.DiscreteGaussian(3.0))
            >>> LWE.dual(params)
            rop: ≈2^116.5, mem: ≈2^73.2, m: 1140, β: 298, d: 2164, ↻: 1, tag: dual
            >>> dual_hybrid(params)
            rop: ≈2^116.2, mem: ≈2^100.4, m: 1137, β: 297, d: 2155, ↻: 1, ζ: 6, tag: dual_hybrid
            >>> dual_hybrid(params, mitm_optimization=True)
            rop: ≈2^160.7, mem: ≈2^156.8, m: 1473, k: 25, ↻: 1, β: 456, d: 2472, ζ: 25, tag: dual_mitm_hybrid

            >>> dual_hybrid(schemes.NTRUHPS2048509Enc)
            rop: ≈2^131.7, mem: ≈2^128.5, m: 436, β: 358, d: 906, ↻: 1, ζ: 38, tag: dual_hybrid

            >>> LWE.dual(schemes.CHHS_4096_67)
            rop: ≈2^206.9, mem: ≈2^137.5, m: ≈2^11.8, β: 616, d: 7779, ↻: 1, tag: dual

            >>> dual_hybrid(schemes.Kyber512, red_cost_model=RC.GJ21, fft=True)
            rop: ≈2^149.8, mem: ≈2^92.1, m: 510, t: 76, β: 399, d: 1000, ↻: 1, ζ: 22, tag: dual_hybrid

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
                log_level=None,
                fft=False,
            ):
                h = params.Xs.hamming_weight
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


class MATZOV:
    """
    See [AC:GuoJoh21]_ and [MATZOV22]_.
    """

    C_prog = 1.0 / (1 - 2.0 ** (-0.292))  # p.37
    C_mul = 32**2  # p.37
    C_add = 5 * 32  # guessing based on C_mul

    @classmethod
    def T_fftf(cls, k, p):
        """
        The time complexity of the FFT in dimension `k` with modulus `p`.

        :param k: Dimension
        :param p: Modulus ≥ 2

        """
        return cls.C_mul * k * p ** (k + 1)  # Theorem 7.6, p.38

    @classmethod
    def T_tablef(cls, D):
        """
        Time complexity of updating the table in each iteration.

        :param D: Number of nonzero entries

        """
        return 4 * cls.C_add * D  # Theorem 7.6, p.39

    @classmethod
    def Nf(cls, params, m, beta_bkz, beta_sieve, k_enum, k_fft, p):
        """
        Required number of samples to distinguish with advantage.

        :param params: LWE parameters
        :param m:
        :param beta_bkz: Block size used for BKZ reduction
        :param beta_sieve: Block size used for sampling
        :param k_enum: Guessing dimension
        :param k_fft: FFT dimension
        :param p: FFT modulus

        """
        mu = 0.5
        k_lat = params.n - k_fft - k_enum  # p.15

        # p.39
        lsigma_s = (
            params.Xe.stddev ** (m / (m + k_lat))
            * (params.Xs.stddev * params.q) ** (k_lat / (m + k_lat))
            * sqrt(4 / 3.0)
            * sqrt(beta_sieve / 2 / pi / e)
            * deltaf(beta_bkz) ** (m + k_lat - beta_sieve)
        )

        # p.29, we're ignoring O()
        N = (
            exp(4 * (lsigma_s * pi / params.q) ** 2)
            * exp(k_fft / 3.0 * (params.Xs.stddev * pi / p) ** 2)
            * (k_enum * cls.Hf(params.Xs) + k_fft * log(p) + log(1 / mu))
        )

        return RR(N)

    @staticmethod
    def Hf(Xs):
        return RR(
            1 / 2 + log(sqrt(2 * pi) * Xs.stddev) + log(coth(pi**2 * Xs.stddev**2))
        ) / log(2.0)

    @classmethod
    def cost(
        cls,
        beta,
        params,
        m=None,
        p=2,
        k_enum=0,
        k_fft=0,
        beta_sieve=None,
        red_cost_model=red_cost_model_default,
    ):
        """
        Theorem 7.6

        """

        if m is None:
            m = params.n

        k_lat = params.n - k_fft - k_enum  # p.15

        # We assume here that β_sieve ≈ β
        N = cls.Nf(
            params,
            m,
            beta,
            beta_sieve if beta_sieve else beta,
            k_enum,
            k_fft,
            p,
        )

        rho, T_sample, _, beta_sieve = red_cost_model.short_vectors(
            beta, N=N, d=k_lat + m, sieve_dim=beta_sieve
        )

        H = cls.Hf(params.Xs)

        coeff = 1 / (1 - exp(-1 / 2 / params.Xs.stddev**2))
        tmp_alpha = pi**2 * params.Xs.stddev**2
        tmp_a = exp(8 * tmp_alpha * exp(-2 * tmp_alpha) * tanh(tmp_alpha)).n(30)
        T_guess = coeff * (
            ((2 * tmp_a / sqrt(e)) ** k_enum)
            * (2 ** (k_enum * H))
            * (cls.T_fftf(k_fft, p) + cls.T_tablef(N))
        )

        cost = Cost(rop=T_sample + T_guess, problem=params)
        cost["red"] = T_sample
        cost["guess"] = T_guess
        cost["beta"] = beta
        cost["p"] = p
        cost["zeta"] = k_enum
        cost["t"] = k_fft
        cost["beta_"] = beta_sieve
        cost["N"] = N
        cost["m"] = m

        cost.register_impermanent({"β'": False, "ζ": False, "t": False}, rop=True, p=False, N=False)
        return cost

    def __call__(
        self,
        params: LWEParameters,
        red_cost_model=red_cost_model_default,
        log_level=1,
    ):
        """
        Optimizes cost of dual attack as presented in [MATZOV22]_.

        See also [AC:GuoJoh21]_.

        :param params: LWE parameters
        :param red_cost_model: How to cost lattice reduction

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``red``: Number of word operations in lattice reduction and
                   short vector sampling.
        - ``guess``: Number of word operations in guessing and FFT.
        - ``β``: BKZ block size.
        - ``ζ``: Number of guessed coordinates.
        - ``t``: Number of coordinates in FFT part mod `p`.
        - ``d``: Lattice dimension.

        """
        params = params.normalize()

        for p in early_abort_range(2, params.q):
            for k_enum in early_abort_range(0, params.n, 10):
                for k_fft in early_abort_range(0, params.n - k_enum[0], 10):
                    # RC.ADPS16(1754, 1754) ~ 2^(512)
                    with local_minimum(40, min(params.n, 1754), log_level=log_level + 4) as it:
                        for beta in it:
                            cost = self.cost(
                                beta,
                                params,
                                p=p[0],
                                k_enum=k_enum[0],
                                k_fft=k_fft[0],
                                red_cost_model=red_cost_model,
                            )
                            it.update(cost)
                        Logging.log(
                            "dual",
                            log_level + 3,
                            f"t: {k_fft[0]}, {repr(it.y)}",
                        )
                        k_fft[1].update(it.y)
                Logging.log("dual", log_level + 2, f"ζ: {k_enum[0]}, {repr(k_fft[1].y)}")
                k_enum[1].update(k_fft[1].y)
            Logging.log("dual", log_level + 1, f"p:{p[0]}, {repr(k_enum[1].y)}")
            p[1].update(k_enum[1].y)
            # if t == 0 then p is irrelevant, so we early abort that loop if that's the case once we hit t==0 twice.
            if p[1].y["t"] == 0 and p[0] > 2:
                break
        Logging.log("dual", log_level, f"{repr(p[1].y)}")
        return p[1].y


matzov = MATZOV()


def dual(
    params: LWEParameters,
    success_probability: float = 0.99,
    red_cost_model=red_cost_model_default,
):
    """
    Dual attack as in [PQCBook:MicReg09]_.

    :param params: LWE parameters.
    :param success_probability: The success probability to target.
    :param red_cost_model: How to cost lattice reduction.

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
    mitm_optimization=False,
    opt_step=8,
    fft=False,
):
    """
    Dual hybrid attack from [INDOCRYPT:EspJouKha20]_.

    :param params: LWE parameters.
    :param success_probability: The success probability to target.
    :param red_cost_model: How to cost lattice reduction.
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
        opt_step=opt_step,
        fft=fft,
    )
    if mitm_optimization:
        ret["tag"] = "dual_mitm_hybrid"
    else:
        ret["tag"] = "dual_hybrid"
    return ret
