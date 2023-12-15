# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sage.all import log, ceil


@dataclass
class SISParameters:
    """The parameters for a Short Integer Solution problem instance."""

    n: int  #: the length of the SIS output.
    q: int  #: the modulus of the space Z/qZ of integers the LWE samples are in.
    length_bound: float  #: The length (in the norm specified below) of an admissable solution.

    #: the length of the SIS preimage.
    #: optionally `sage.all.oo` for allowing infinitely many samples.
    m: int = None
    norm: int = 2  #: The norm to use for measuring length (l-p norm) use oo for infinity norm.

    tag: str = None  #: a name for the patameter set

    def __post_init__(self, **kwds):
        if not self.m:
            self.m = ceil(self.n*log(self.q))  #: Set m to be the minimum required for a solution to exist.

    @property
    def _homogeneous(self):
        return False

    def updated(self, **kwds):  # TODO Add docstrings for SIS scheme parameters based on Dilithium.
        """
        Return a new set of parameters updated according to ``kwds``.

        :param kwds: We set ``key`` to ``value`` in the new set of parameters.

        EXAMPLE::

            >>> from estimator import *

        """
        d = dict(self.__dict__)
        d.update(kwds)
        return SISParameters(**d)