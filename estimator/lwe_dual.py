from functools import partial
from dataclasses import replace

from sage.all import oo, ceil, sqrt, log, cached_function, exp
from .reduction import delta as deltaf
from .reduction import cost as costf
from .reduction import ADPS16, BDGL16
from .reduction import LLL
from .util import binary_search, robust_bin_search
from .cost import Cost
from .lwe import LWEParameters
from .prob import drop as prob_drop
from .prob import amplify as prob_amplify
from .io import Logging
from .conf import red_cost_model_default, default_mitm_opt
from .errors import OutOfBoundsError
from .nd import NoiseDistribution
from .lwe_brute_force import exhaustive_search, mitm


class DualHybrid:

    full_sieves = [ADPS16, BDGL16]

    def __init__(self, opt_step=1):
        self.opt_step = opt_step

    @staticmethod
    @cached_function
    def dual_reduce(
        delta_0: float,
        params: LWEParameters,
        zeta: int = 0,
        h1: int = 0,
        scaling_factor: float = 1.,
        log_level=None,
    ):
        """
        Produce new LWE sample using a dual vector on first n-zeta coordinates of the secret. The
        length of the dual vector is given by delta_0 in root Hermite form and using a possible
        scaling factor, i.e. |v| = scaling_factor * delta_0^d * q^((n-zeta)/d).

        :param delta_0: Length of the vector in root Hermite form
        :param params: LWE parameters
        :param zeta: Dimension ζ ≥ 0 of new LWE instance
        :param h1: Number of non-zero components of the secret of the new LWE instance
        :param scaling_factor: Factor introduced by obtaining multiple dual vectors

        .. note :: This function assumes that the instance is normalized.
        """
        if not 0 <= zeta <= params.n:
            raise OutOfBoundsError(f"Splitting dimension {zeta} must be between 0 and n={params.n}.")

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

        # see if we have optimally many samples (as in [ia.cr/2020/515]) available
        m_ = max(1, ceil(sqrt(red_Xs.n * log(c) / log(delta_0))) - red_Xs.n)
        m_ = min(params.m, m_)

        # Compute new noise as in [ia.cr/2020/515]
        sigma_ = scaling_factor * red_Xs.stddev * delta_0 ** (m_ + red_Xs.n) / c ** (m_ / (m_ + red_Xs.n))
        slv_Xe = NoiseDistribution.DiscreteGaussian(params.q * sigma_)

        slv_params = LWEParameters(
                        n=zeta,
                        q=params.q,
                        Xs=slv_Xs,
                        Xe=slv_Xe,
                    )
        return slv_params, m_

    @staticmethod
    @cached_function
    def cost(
        solver,
        params: LWEParameters,
        beta: int,
        zeta: int = 0,
        h1: int = 0,
        success_probability: float = .99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        log_level=None,
    ):
        """
        Computes the cost of the dual hybrid attack that dual reduces the LWE instance and then
        uses the given solver to solve the reduced instance.

        :param solver: Algorithm for solving the reduced instance
        :param params: LWE parameters
        :param beta: Blocksize used to produce short dual vectors for dual reduction
        :param zeta: Dimension ζ ≥ 0 of new LWE instance
        :param h1: Number of non-zero components of the secret of the new LWE instance
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction
        :param use_lll: Use LLL calls to produce more small vectors

        .. note :: This function assumes that the instance is normalized. It runs no optimization,
            it merely reports costs.
        """
        Logging.log("dual", log_level, f"costing with beta={beta}, zeta={zeta}, h1={h1}")

        delta_0 = deltaf(beta)

        if red_cost_model in DualHybrid.full_sieves:
            scaling_factor = 4. / 3
        elif use_lll:
            scaling_factor = 2.
        else:
            scaling_factor = 1.

        params_slv, m_ = DualHybrid.dual_reduce(delta_0, params, zeta, h1, scaling_factor, log_level=log_level+1)
        Logging.log("dual", log_level + 1, f"reduced LWE instance: {repr(params_slv)}")

        cost_slv = solver(params_slv, success_probability)
        Logging.log("dual", log_level + 2, f"cost of solver: {repr(cost_slv)}")

        d = m_ + params.n - zeta
        cost_red = costf(red_cost_model, beta, d)
        if red_cost_model in DualHybrid.full_sieves:
            # if we use full sieving, we get many short vectors
            # we compute in logs to avoid overflows in case m
            # or beta happen to be large
            try:
                log_rep = max(0, log(cost_slv["m"]) - (beta / 2) * log(4 / 3))
                cost_red = cost_red.repeat(ceil(exp(log_rep)))
            except OverflowError:
                # if we still get an overflow, m must be huge so
                # can probably be approximated with oo for our
                # purposes
                cost_slv["rop"] = oo
        elif use_lll:
            cost_red["rop"] += cost_slv["m"] * LLL(d, log(params.q, 2))
            cost_red["repetitions"] = cost_slv["m"]
        else:
            cost_red = cost_red.repeat(cost_slv["m"])

        Logging.log("dual", log_level + 2, f"cost of reduction: {repr(cost_red)}")

        total_cost = cost_slv.combine(cost_red)
        total_cost["m"] = m_
        total_cost["rop"] = cost_red["rop"] + cost_slv["rop"]
        total_cost["mem"] = cost_slv["mem"]

        if d < params.n - zeta:
            raise RuntimeError(f"{d} < {params.n - zeta}, {params.n}, {zeta}, {m_}")
        total_cost["d"] = d

        Logging.log("dual", log_level, f"total cost: {repr(total_cost)}")

        rep = 1
        if params.Xs.is_sparse:
            h = params.Xs.get_hamming_weight(params.n)
            probability = prob_drop(params.n, h, zeta, h1)
            rep = prob_amplify(success_probability, probability)
        # don't need more samples to re-run attack, since we may
        # just guess different components of the secret
        return total_cost.repeat(times=rep, select={"m" : False})

    def optimize_blocksize(
        self,
        solver,
        params: LWEParameters,
        zeta: int = 0,
        h1: int = 0,
        success_probability: float = .99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        log_level=None,
    ):
        """
        Optimizes the cost of the dual hybrid attack over the blocksize beta.

        :param solver: Algorithm for solving the reduced instance
        :param params: LWE parameters
        :param zeta: Dimension ζ ≥ 0 of new LWE instance
        :param h1: Number of non-zero components of the secret of the new LWE instance
        :param success_probability: The success probability to target
        :param red_cost_model: How to cost lattice reduction
        :param use_lll: Use LLL calls to produce more small vectors

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
            log_level=log_level)
        # don't have a reliable upper bound for beta
        # we choose n - k arbitrarily and adjust later if
        # necessary
        beta_upper = max(params.n - zeta, 40)
        beta = beta_upper
        while beta == beta_upper:
            beta_upper *= 2
            cost = robust_bin_search(f, 2, beta_upper, "beta", step=self.opt_step)
            beta = cost["beta"]

        cost["zeta"] = zeta
        if params.Xs.is_sparse:
            cost["h1"] = h1
        return cost

    def __call__(
        self,
        solver,
        params: LWEParameters,
        success_probability: float = .99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        opt_step=2,
        log_level=1,
    ):
        """
        Optimizes the cost of the dual hybrid attack (using the given solver) over
        all attack parameters: blocksize beta, splitting dimension ζ, and
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
        :param use_lll: use LLL calls to produce more small vectors
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
        - When ζ > 1 0 and ``solver`` is ``exhaustive_search`` this function estimates
            the hybrid attack as given in [INDOCRYPT:EspJouKha20]_
        - When ζ > 1 0 and ``solver`` is ``mitm`` this function estimates the dual MITM
            hybrid attack roughly following [https://ia.cr/2019/1114]_

        EXAMPLE::

            >>> from estimator import *
            >>> params = LWEParameters(n=1024, q = 2**32, Xs=ND.Uniform(0,1), Xe=ND.DiscreteGaussian(3.0))
            >>> dual(params)
            rop: ≈2^115.5, mem: ≈2^70.0, m: 1018, red: ≈2^115.4, δ: 1.005021, β: 284, d: 2041, ↻: ≈2^69.0
            >>> dual_hybrid(params)
            rop: ≈2^111.3, mem: ≈2^106.4, m: 983, red: ≈2^111.2, δ: 1.005204, β: 269, d: 1957, ↻: ≈2^56.4, ζ: 50
            >>> dual_mitm_hybrid(params)
            rop: ≈2^141.1, mem: ≈2^139.1, m: 1189, k: 132, ↻: 139, red: ≈2^140.8, δ: 1.004164, β: 375, d: 2021, ζ: 192
            >>> dual_mitm_hybrid(params, mitm_optimization="numerical")
            rop: ≈2^140.6, m: 1191, k: 128, mem: ≈2^136.0, ↻: 133, red: ≈2^140.2, δ: 1.004179, β: 373, d: 11, ζ: 163

            >>> from dataclasses import replace
            >>> params = replace(params, Xs=ND.SparseTernary(params.n, 32))
            >>> dual(params)
            rop: ≈2^112.8, mem: ≈2^64.0, m: 953, red: ≈2^112.7, δ: 1.005178, β: 271, d: 1976, ↻: ≈2^65.0
            >>> dual_hybrid(params)
            rop: ≈2^97.8, mem: ≈2^81.9, m: 730, red: ≈2^97.4, δ: 1.006813, β: 175, d: 1453, ↻: ≈2^36.3, ζ: 301, h1: 8
            >>> dual_mitm_hybrid(params)
            rop: ≈2^103.4, mem: ≈2^81.5, m: 724, k: 310, ↻: ≈2^27.3, red: ≈2^102.7, δ: 1.006655, β: 182, ...
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
        self.opt_step = opt_step
        Logging.log("dual", log_level, f"costing LWE instance: {repr(params)}")

        params = params.normalize()

        if params.Xs.is_sparse:
            Cost.register_impermanent(h1=False)

            def _optimize_blocksize(
                solver,
                params: LWEParameters,
                zeta: int = 0,
                success_probability: float = .99,
                red_cost_model=red_cost_model_default,
                use_lll=True,
                log_level=None
            ):
                f = partial(
                    self.optimize_blocksize,
                    solver=solver,
                    params=params,
                    zeta=zeta,
                    success_probability=success_probability,
                    red_cost_model=red_cost_model,
                    use_lll=use_lll,
                    log_level=log_level + 1
                    )
                h = params.Xs.get_hamming_weight(params.n)
                h1_min = max(0, h - (params.n - zeta))
                h1_max = min(zeta, h) - 1       # subtracting 1 as workaround for issue #6
                Logging.log("dual", log_level, f"Optimizing over h1 in [{h1_min},{h1_max}] (zeta={zeta})")
                return binary_search(f, h1_min, h1_max, "h1")
        else:
            _optimize_blocksize = self.optimize_blocksize

        f = partial(
            _optimize_blocksize,
            solver=solver,
            params=params,
            success_probability=success_probability,
            red_cost_model=red_cost_model,
            use_lll=use_lll,
            log_level=log_level+1,
            )

        return robust_bin_search(f, 1, params.n - 1, "zeta", step=self.opt_step)


DH = DualHybrid(opt_step=2)


def dual(
        params: LWEParameters,
        success_probability: float = .99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
):
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
        solver=exhaustive_search,
        params=params,
        zeta=1,
        h1=0,
        success_probability=success_probability,
        red_cost_model=red_cost_model,
        use_lll=use_lll,
        log_level=1
        )
    del ret["zeta"]
    if hasattr(ret, "h1"):
        del ret["h1"]
    return ret


def dual_hybrid(
        params: LWEParameters,
        success_probability: float = .99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        opt_step=2
):
    return DH(
        solver=exhaustive_search,
        params=params,
        success_probability=success_probability,
        red_cost_model=red_cost_model,
        use_lll=use_lll,
        opt_step=opt_step
        )


def dual_mitm_hybrid(
        params: LWEParameters,
        success_probability: float = .99,
        red_cost_model=red_cost_model_default,
        use_lll=True,
        mitm_optimization=default_mitm_opt,
        opt_step=2
):
    solver = partial(mitm, optimization=mitm_optimization)
    return DH(
        solver=solver,
        params=params,
        success_probability=success_probability,
        red_cost_model=red_cost_model,
        use_lll=use_lll,
        opt_step=opt_step
        )
