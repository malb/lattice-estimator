# -*- coding: utf-8 -*-
from dataclasses import dataclass

from .conf import ntru_fatigue_lb
from .errors import InsufficientSamplesError
from .lwe_parameters import LWEParameters


@dataclass
class NTRUParameters(LWEParameters):
    """The parameters for an NTRU problem instance. The estimator treats regular NTRU parameters as similar
    to LWE, but requires different estimation methodology for overstrethed parameters.

    :param ntru_type:
        Specifies the type of NTRU instance the parameters represent. Currently supported
        types are, "matrix" for general matrix NTRU, "circulant" for circulant NTRU, "fixed" for circulant
        NTRU with a fixed geometry.

    """

    ntru_type: str = "matrix"

    def __post_init__(self, **kwds):
        super().__post_init__()
        # set m = n
        self.m = self.n

    @property
    def possibly_overstretched(self):
        if self.q >= ntru_fatigue_lb(self.n):
            return True

        return False

    @property
    def homogeneous(self):
        return True

    def normalize(self):
        """
        EXAMPLES:

        We perform the normal form transformation if χ_e < χ_s and we got the samples::
        For NTRU, m = n so we swap the secret and the noise::

            >>> from estimator import *
            >>> Xs=ND.DiscreteGaussian(2.0)
            >>> Xe=ND.DiscreteGaussian(1.58)
            >>> NTRU.Parameters(n=512, q=8192, Xs=Xs, Xe=Xe, m=512).normalize()
            NTRUParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=2.00), m=512, tag=None, ntru_type='matrix')

        """
        if self.m < 1:
            raise InsufficientSamplesError(f"m={self.m} < 1")

        # swap secret and noise
        # TODO: this is somewhat arbitrary
        if self.Xe < self.Xs and self.m < 2 * self.n:
            return NTRUParameters(n=self.n, q=self.q, Xs=self.Xe, Xe=self.Xs, m=self.n,
                                  tag=self.tag, ntru_type=self.ntru_type)

        # nothing to do
        return self

    def updated(self, **kwds):
        """
        Return a new set of parameters updated according to ``kwds``.

        :param kwds: We set ``key`` to ``value`` in the new set of parameters.

        EXAMPLE::

            >>> from estimator import *
            >>> schemes.NTRUHPS2048509Enc  #doctest: +ELLIPSIS
            NTRUParameters(n=508, q=2048, Xs=D(σ=0.82), Xe=D(σ=0.71), m=508, tag='NTRUHPS2048509Enc', ntru_type='ma...
            >>> schemes.NTRUHPS2048509Enc.possibly_overstretched
            False

            >>> schemes.NTRUHPS2048509Enc.updated(q=16536)  #doctest: +ELLIPSIS
            NTRUParameters(n=508, q=16536, Xs=D(σ=0.82), Xe=D(σ=0.71), ..., tag='NTRUHPS2048509Enc', ntru_type='matrix')
            >>> schemes.NTRUHPS2048509Enc.updated(q=16536).possibly_overstretched
            True
        """
        d = dict(self.__dict__)
        d.update(kwds)
        return NTRUParameters(**d)

    def amplify_m(self, m):
        """
        Return an NTRU instance parameters with ``m`` samples produced from the samples in this instance.

        :param m: New number of samples.

        """
        raise NotImplementedError("Rerandomizing NTRU instances is not supported yet.")

    def switch_modulus(self):
        """
        Apply modulus switching and return new instance.

        See [JMC:AlbPlaSco15]_ for details.
        """
        raise NotImplementedError("Modulus Switching for NTRU not supported yet.")

    def __hash__(self):
        return hash((self.n, self.q, self.Xs, self.Xe, self.m, self.tag, self.ntru_type))
