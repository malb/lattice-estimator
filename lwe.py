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

    def __hash__(self):
        return hash((self.n, self.q, self.Xs, self.Xe, self.m, self.tag))
