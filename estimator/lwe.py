# -*- coding: utf-8 -*-
from sage.all import oo
from dataclasses import dataclass
from .nd import NoiseDistribution
from .errors import InsufficientSamplesError


@dataclass
class LWEParameters:
    n: int
    q: int
    Xs: NoiseDistribution
    Xe: NoiseDistribution
    m: int = oo
    tag: str = None

    def __post_init__(self, **kwds):
        self.Xs.n = self.n
        if self.m < oo:
            self.Xe.n = self.m

    def normalize(self):
        """
        EXAMPLES:

        We perform the normal form transformation if χ_e < χ_s and we got the samples::

            >>> from estimator import *
            >>> Xs=ND.DiscreteGaussian(2.0)
            >>> Xe=ND.DiscreteGaussian(1.58)
            >>> LWEParameters(n=512, q=8192, Xs=Xs, Xe=Xe).normalize()
            LWEParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=1.58), m=+Infinity, tag=None)

        If m = n, we swap the secret and the noise::

            >>> from estimator import *
            >>> Xs=ND.DiscreteGaussian(2.0)
            >>> Xe=ND.DiscreteGaussian(1.58)
            >>> LWEParameters(n=512, q=8192, Xs=Xs, Xe=Xe, m=512).normalize()
            LWEParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=2.00), m=512, tag=None)

        """
        if self.m < 1:
            raise InsufficientSamplesError(f"m={self.m} < 1")

        # Normal form transformation
        if self.Xe < self.Xs and self.m >= 2 * self.n:
            return LWEParameters(
                n=self.n, q=self.q, Xs=self.Xe, Xe=self.Xe, m=self.m - self.n, tag=self.tag
            )
        # swap secret and noise
        if self.Xe < self.Xs and self.m == self.n:
            return LWEParameters(n=self.n, q=self.q, Xs=self.Xe, Xe=self.Xs, m=self.n, tag=self.tag)

        # nothing to do
        return self

    def updated(self, **kwds):
        """
        Return a new set of parameters updated according to ``kwds``.

        :param kwds: We set ``key`` to ``value`` in the new set of parameters.

        EXAMPLE::

            >>> from estimator import *
            >>> Kyber512
            LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.00), m=1024, tag='Kyber 512')
            >>> Kyber512.updated(m=1337)
            LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.00), m=1337, tag='Kyber 512')

        """
        d = dict(self.__dict__)
        d.update(kwds)
        return LWEParameters(**d)

    def __hash__(self):
        return hash((self.n, self.q, self.Xs, self.Xe, self.m, self.tag))
