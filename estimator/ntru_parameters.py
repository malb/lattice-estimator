# -*- coding: utf-8 -*-
from dataclasses import dataclass
from copy import copy
from sys import stderr

from sage.all import oo, binomial, log, sqrt, ceil, exp

from .conf import ntru_fatigue_lb, ntru_fatigue_ub
from .nd import NoiseDistribution
from .errors import InsufficientSamplesError
from .lwe_parameters import LWEParameters
from .overstretched_ntru import find_fatigue


@dataclass
class NTRUParameters(LWEParameters):
    """The parameters for an NTRU problem instance. The estimator treats regular NTRU parameters as similar
    to LWE, but requires different estimation methodology for overstrethed parameters. """

    ntru_type: str = "matrix"

    def __post_init__(self, **kwds):
        super().__post_init__()
        
        # Use lower bound on fatigue point to inform user on possible overstretched parameter set
        # fatigue_point = find_fatigue(self.n, sk_variance=self.Xs.stddev**2, ntru=self.ntru_type)
        if self.q >= ntru_fatigue_lb(self.n):
            self.possibly_overstretched = True

            # TODO: Make below print statement use proper logging procedures.
            print(f"NOTE: NTRU with n = {self.n}, q = {self.q} is potentially overstretched. Primal attack estimation will include fatigue point estimation.", file=stderr)

        else:
            self.possibly_overstretched = False


# TODO: Below are LWE specific transformations. Which ones apply to ntru?

    def normalize(self):
        """
        EXAMPLES:

        We perform the normal form transformation if χ_e < χ_s and we got the samples::
        For NTRU, m = n so we swap the secret and the noise::

            >>> from estimator import *
            >>> Xs=ND.DiscreteGaussian(2.0)
            >>> Xe=ND.DiscreteGaussian(1.58)
            >>> NTRU.Parameters(n=512, q=8192, Xs=Xs, Xe=Xe, m=512).normalize()
            NTRUParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=2.00), m=512, tag=None, ntru_type="matrix")

        """
        if self.m < 1:
            raise InsufficientSamplesError(f"m={self.m} < 1")

        # swap secret and noise
        # TODO: this is somewhat arbitrary
        if self.Xe < self.Xs and self.m < 2 * self.n:
            return NTRUParameters(n=self.n, q=self.q, Xs=self.Xe, Xe=self.Xs, m=self.n, tag=self.tag, ntru_type=self.ntru_type)

        # nothing to do
        return self

    def updated(self, **kwds):
        """
        Return a new set of parameters updated according to ``kwds``.

        :param kwds: We set ``key`` to ``value`` in the new set of parameters.

        EXAMPLE::

            >>> from estimator import *
            >>> schemes.Kyber512
            LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.22), m=512, tag='Kyber 512')
            >>> schemes.Kyber512.updated(m=1337)
            LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.22), m=1337, tag='Kyber 512')

        """
        d = dict(self.__dict__)
        d.update(kwds)
        return NTRUParameters(**d)

    def amplify_m(self, m):
        """
        Return a LWE instance parameters with ``m`` samples produced from the samples in this instance.

        :param m: New number of samples.

        EXAMPLE::

            >>> from sage.all import binomial, log
            >>> from estimator import *
            >>> schemes.Kyber512
            LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.22), m=512, tag='Kyber 512')
            >>> schemes.Kyber512.amplify_m(2**100)
            LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=4.58), m=..., tag='Kyber 512')

        We can produce 2^100 samples by random ± linear combinations of 12 vectors::

            >>> float(sqrt(12.)), float(log(binomial(1024, 12) , 2.0)) + 12
            (3.46..., 103.07...)

        """
        raise NotImplementedError("Rerandomizing NTRU instances is not supported yet.")

    def switch_modulus(self):
        """
        Apply modulus switching and return new instance.

        See [JMC:AlbPlaSco15]_ for details.

        EXAMPLE::

            >>> from estimator import *
            >>> LWE.Parameters(n=128, q=7681, Xs=ND.UniformMod(3), Xe=ND.UniformMod(11)).switch_modulus()
            LWEParameters(n=128, q=5289, Xs=D(σ=0.82), Xe=D(σ=3.08), m=+Infinity, tag=None)

        """
        # n = self.Xs.density * len(self.Xs)

        # # n uniform in -(0.5,0.5) ± stddev(χ_s)
        # Xr_stddev = sqrt(n / 12) * self.Xs.stddev  # rounding noise
        # # χ_r == p/q ⋅ χ_e # we want the rounding noise match the scaled noise
        # p = ceil(Xr_stddev * self.q / self.Xe.stddev)

        # scale = float(p) / self.q

        # # there is no point in scaling if the improvement is eaten up by rounding noise
        # if scale > 1 / sqrt(2):
        #     return self

        # return LWEParameters(
        #     self.n,
        #     p,
        #     Xs=self.Xs,
        #     Xe=NoiseDistribution.DiscreteGaussian(sqrt(2) * self.Xe.stddev * scale),
        #     m=self.m,
        #     tag=f"{self.tag},scaled" if self.tag else None,
        # )
        raise NotImplementedError("Modulus Switching for NTRU not supported yet.")

    def __hash__(self):
        return hash((self.n, self.q, self.Xs, self.Xe, self.m, self.tag, self.ntru_type))
