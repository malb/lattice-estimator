from functools import partial

from sage.all import sqrt, pi, log, exp, RR, ZZ, oo, round, e, binomial, ceil, floor
from .cost import Cost
from .lwe import LWEParameters
from .errors import InsufficientSamplesError
from .nd import sigmaf
from .conf import default_mitm_opt
from .util import binary_search, log2
from .prob import amplify


class ExhaustiveSearch:

    def __call__(
        self,
        params: LWEParameters,
        success_probability=0.99,
        quantum: bool = False
    ):
        """
        Estimate cost of solving LWE via exhaustive search.

        :param params: LWE parameters
        :param success_probability: the targeted success probability
        :param quantum: use estimate for quantum computer (we simply take the square root of the search space)
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: memory requirement in integers mod q.
        - ``m``: Required number of samples to distinguish the correct solution with high probability.

        EXAMPLE::

            >>> from estimator import *
            >>> params = LWEParameters(n=64, q=2**40, Xs=ND.UniformMod(2), Xe=ND.DiscreteGaussian(3.2))
            >>> exhaustive_search(params)
            rop: ≈2^73.6, mem: ≈2^72.6, m: 397.198
            >>> params = LWEParameters(n=1024, q=2**40, Xs=ND.SparseTernary(n=1024, p=32), Xe=ND.DiscreteGaussian(3.2))
            >>> exhaustive_search(params)
            rop: ≈2^417.3, mem: ≈2^416.3, m: ≈2^11.2

        """
        params = LWEParameters.normalize(params)

        # there are two stages: enumeration and distinguishing, so we split up the success_probability
        probability = sqrt(success_probability)

        size = params.Xs.support_size(n=params.n, fraction=probability)

        if quantum:
            size = size.sqrt()

        # set m according to [ia.cr/2020/515]
        sigma = params.Xe.stddev / params.q
        m_required = RR(8 * exp(4 * pi * pi * sigma * sigma) * (log(size) - log(log(1 / probability))))

        if params.m < m_required:
            raise InsufficientSamplesError(
                f"Exhaustive search: Need {m_required} samples but only {params.m} available."
                )
        else:
            m = m_required

        # we can compute A*s for all candidate s in time 2*size*m using
        # (the generalization [ia.cr/2021/152] of) the recursive algorithm
        # from [ia.cr/2020/515]
        cost = 2 * size * m

        ret = Cost(rop=cost, mem=cost / 2, m=m)
        return ret


exhaustive_search = ExhaustiveSearch()


class MITM:
    t = 2       # cut-off for discret Gaussian distributions
    p = 1 - 2 * exp(-4 * pi)    # probability that a coefficient falls within the cut-off

    locality = .05

    def X_range(self, nd):
        a, b = nd.bounds
        p = 1.
        if a == -oo or b == oo:
            # Looks like nd is Gaussian
            # -> we'll treat it as bounded (with probability self.p)
            b = 2 * self.t * sigmaf(nd.stddev) + 1
            a = -b
            p = self.p
        return b - a + 1, p

    def local_range(self, center):
        return ZZ(floor((1 - self.locality) * center)), ZZ(ceil((1 + self.locality) * center))

    def mitm_analytical(
        self,
        params: LWEParameters,
        success_probability=0.99
    ):
        nd_rng, nd_p = self.X_range(params.Xe)
        delta = nd_rng / params.q        # possible error range scaled

        sd_rng, sd_p = self.X_range(params.Xs)

        # determine the number of elements in the tables depending on splitting dim
        n = params.n
        k = round(n / (2 - delta))
        # we could now call self.cost with this k, but using our model below seems
        # about 3x faster and reasonably accurate

        if params.Xs.is_sparse:
            h = params.Xs.get_hamming_weight(n=params.n)
            split_h = round(h * k / n)
            success_probability_ = binomial(k, split_h) * binomial(n - k, h - split_h) / binomial(n, h)

            logT = RR(h * (log2(n) - log2(h) + log2(sd_rng - 1) + log2(e))) / (2 - delta)
            logT -= RR(log2(h) / 2)
            logT -= RR(h * h * log2(e) / (2 * n * (2 - delta) ** 2))
        else:
            success_probability_ = 1.
            logT = k * log(sd_rng, 2)

        m_ = max(1, round(logT + log(logT, 2)))
        if params.m < m_:
            raise InsufficientSamplesError(
                f"MITM: Need {m_} samples but only {params.m} available."
                )

        # since m = logT + loglogT and rop = T*m, we have rop=2^m
        ret = Cost(rop=2 ** m_, mem=2 ** logT * m_, m=m_, k=ZZ(k))
        repeat = amplify(success_probability, sd_p ** n * nd_p ** m_ * success_probability_)
        return ret.repeat(times=repeat)

    def cost(
        self,
        params: LWEParameters,
        k:  int,
        success_probability=0.99,
    ):
        nd_rng, nd_p = self.X_range(params.Xe)
        delta = nd_rng / params.q        # possible error range scaled

        sd_rng, sd_p = self.X_range(params.Xs)
        n = params.n

        if params.Xs.is_sparse:
            h = params.Xs.get_hamming_weight(n=n)

            # we assume the hamming weight to be distributed evenly across the two parts
            # if not we can rerandomize on the coordinates and try again -> repeat
            split_h = round(h * k / n)
            size_tab = RR((sd_rng - 1) ** split_h * binomial(k, split_h))
            size_sea = RR((sd_rng - 1) ** (h - split_h) * binomial(n - k, h - split_h))
            success_probability_ = binomial(k, split_h) * binomial(n - k, h - split_h) / binomial(n, h)
        else:
            size_tab = sd_rng ** k
            size_sea = sd_rng ** (n - k)
            success_probability_ = 1

        # for search we effectively build a second table and for each entry, we expect
        # 2^( m * 4 * B / q) = 2^(delta * m) table look ups + a l_oo computation (costing m)
        # for every hit in the table (which has probability T/2^m)
        def cost_search(m):
            return (m, size_sea * (2 * m + 2 ** (delta * m) * (1 + size_tab * m / 2 ** m)))

        # we set m such that it approximately minimizes the search cost per query as
        # a reasonable starting point and then optimize around it
        m_ = ceil(max(log2(size_tab) + log2(log2(size_tab)), 1))
        a, b = self.local_range(m_)
        m_, cost_search = binary_search(
                                cost_search,
                                a,
                                b,
                                "m",
                                smallerf=lambda x, best: x[1] <= best[1]
                                )
        m_ = min(m_, params.m)

        # building the table costs 2*T*m using the generalization [ia.cr/2021/152] of
        # the recursive algorithm from [ia.cr/2020/515]
        cost_table = size_tab * 2 * m_

        ret = Cost(rop=(cost_table + cost_search), m=m_, k=k)
        ret["mem"] = (size_tab * (k + m_) + size_sea * (n - k + m_))
        repeat = amplify(success_probability, sd_p ** n * nd_p ** m_ * success_probability_)
        return ret.repeat(times=repeat)

    def __call__(
        self,
        params: LWEParameters,
        success_probability=0.99,
        optimization=default_mitm_opt
    ):
        """
        Estimate cost of solving LWE via Meet-In-The-Middle attack.

        :param params: LWE parameters
        :param success_probability: the targeted success probability
        :param model: Either "analytical" (faster, default) or "numerical" (more accurate)
        :return: A cost dictionary

        The returned cost dictionary has the following entries:

        - ``rop``: Total number of word operations (≈ CPU cycles).
        - ``mem``: memory requirement in integers mod q.
        - ``m``: Required number of samples to distinguish the correct solution with high probability.
        - ``k``: Splitting dimension.
        - ``↻``: Repetitions required to achieve targeted success probability

        EXAMPLE::

            >>> from estimator import *
            >>> params = LWEParameters(n=64, q=2**40, Xs=ND.UniformMod(2), Xe=ND.DiscreteGaussian(3.2))
            >>> mitm(params)
            rop: ≈2^37.0, mem: ≈2^37.2, m: 37, k: 32, ↻: 1
            >>> mitm(params, optimization="numerical")
            rop: ≈2^39.2, m: 36, k: 32, mem: ≈2^39.1, ↻: 1
            >>> params = LWEParameters(n=1024, q=2**40, Xs=ND.SparseTernary(n=1024, p=32), Xe=ND.DiscreteGaussian(3.2))
            >>> mitm(params)
            rop: ≈2^215.4, mem: ≈2^210.2, m: ≈2^13.1, k: 512, ↻: 43
            >>> mitm(params, optimization="numerical")
            rop: ≈2^216.0, m: ≈2^13.1, k: 512, mem: ≈2^211.4, ↻: 43

        """
        Cost.register_impermanent(rop=True, mem=False, m=True, k=False)

        params = LWEParameters.normalize(params)

        nd_rng, _ = self.X_range(params.Xe)
        if nd_rng >= params.q:
            # MITM attacks cannot handle an error this large.
            return Cost(rop=oo, mem=oo, m=0, k=0)

        if "analytical" in optimization:
            return self.mitm_analytical(params=params, success_probability=success_probability)
        elif "numerical" in optimization:
            f = partial(self.cost, params=params, success_probability=success_probability)
            ret = binary_search(f, 1, params.n - 1, "k")

            # if the noise is large, the curve might not be convex, so the above minimum
            # is not correct. Interestingly, in these cases, it seems that k=1 might be smallest
            ret1 = f(k=1)
            return min(ret, ret1)
        else:
            raise ValueError("Unknown optimization method for MITM.")


mitm = MITM()
