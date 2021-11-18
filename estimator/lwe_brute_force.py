from sage.all import sqrt, pi, log, exp, RR
from .cost import Cost
from .lwe import LWEParameters
from .errors import InsufficientSamplesError


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
                "Exhaustive Search: Need more than {} samples to recognize correct solution.".format(params.m)
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

# TODO: MITM
