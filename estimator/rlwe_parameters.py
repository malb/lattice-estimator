# -*- coding: utf-8 -*-
import math
from sage.all import oo
from dataclasses import dataclass
from copy import copy

from .lwe_parameters import LWEParameters
from .nd import NoiseDistribution


@dataclass
class RLWEParameters:
    """The parameters for a Ring Learning With Errors problem instance."""

    N: int  #: The degree of the polynomial modulus x^N + 1
    q: int  #: The modulus of the space Z/qZ the coefficients are in.
    Xs: NoiseDistribution  #: the distribution on Z/qZ from which the LWE secret is drawn.
    Xe: NoiseDistribution  #: the distribution on Z/qZ from which the error term is drawn.

    #: the number of LWE samples allowed to an attacker,
    #: optionally `sage.all.oo` for allowing infinitely many samples.
    m: int = oo

    tag: str = None  #: a name for the patameter set

    def __post_init__(self, **kwds):
        self.Xs = copy(self.Xs)
        self.Xs.n = 1  # RLWE instances are a dot product of two polynomials.

        log_N = math.log2(self.N)
        if not log_N.is_integer():
            raise ValueError(
                "RLWE security estimation only supports power-of-2 polynomial "
                f"modulus degrees, but {self.N} (log2={log_N}) was given to "
                "RLWEParameters.")

        if self.m < oo:
            self.Xe = copy(self.Xe)
            self.Xe.n = self.m

    def as_lwe_params(self) -> LWEParameters:
        """Reduce this instance of RLWE to LWE.

    This operates under the assumption that no current known attacks exploit the
    ring structure of RLWE beyond its LWE structure.

    The reduction starts with an RLWE instance with sample

        a(x) = a_0 + a_1 x + ... + a_{N-1}x^{N-1}

    and secret key

        s(x) = s_0 + s_1 x + ... + s_{N-1}x^{N-1},

    where b(x) = a(x)s(x) + m(x) + e(x) mod (x^N + 1).

    Noting that b_i x^i = sum_{j+k=i mod N} a_j s_k x^{j+k},
    (i.e., values have their sign inverted if j+k > N), the reduction converts
    the instance to the following set of LWE samples:

        (a_0,     -a_{N-1}, -a_{N-2}, ..., -a_1),
        (a_1,      a_0,     -a_{N-1}, ..., -a_2),
        (a_2,      a_1,      a_0,     ..., -a_3),
         ...                          ...
        (a_{N-1},  a_{N-2},  a_{N-3}, ...,  a_0),

    With secret s = (s_0, s_1, ..., s_{N-1}),
    and b = (b_0, b_1, ..., b_{N-1}).

    As such, an attack against LWE with dimension N is also an attack against
    RLWE for polynomials of degree N.
    """
        return LWEParameters(
            n=self.N,
            q=self.q,
            Xs=self.Xs,
            Xe=self.Xe,
            m=self.m,
            tag=f"LWE_reduced_from_{self.tag}",
        )

    def __hash__(self):
        return hash((self.N, self.q, self.Xs, self.Xe, self.m, self.tag))
