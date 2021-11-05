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

        :params **kwds: We set ``key`` to ``value`` in the new set of parameters.

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
