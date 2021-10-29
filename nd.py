# -*- coding: utf-8 -*-

from dataclasses import dataclass

from sage.all import parent, RR, RealField, sqrt, pi


def stddevf(sigma):
    """
    Gaussian width parameter σ → standard deviation

    :param sigma: Gaussian width parameter σ

    EXAMPLE::

        sage: from estimator import stddevf
        sage: stddevf(64.0)
        25.532...

        sage: stddevf(64)
        25.532...

        sage: stddevf(RealField(256)(64)).prec()
        256

    """

    try:
        prec = parent(sigma).prec()
    except AttributeError:
        prec = 0
    if prec > 0:
        FF = parent(sigma)
    else:
        FF = RR
    return FF(sigma) / FF(sqrt(2 * pi))


def sigmaf(stddev):
    """
    standard deviation → Gaussian width parameter σ

    :param sigma: standard deviation

    EXAMPLE::

        sage: from estimator import stddevf, sigmaf
        sage: n = 64.0
        sage: sigmaf(stddevf(n))
        64.000...

        sage: sigmaf(RealField(128)(1.0))
        2.5066282746310005024157652848110452530
        sage: sigmaf(1.0)
        2.50662827463100
        sage: sigmaf(1)
        2.50662827463100
        sage: sigmaf(1r)
        2.50662827463100

    """
    RR = parent(stddev)
    #  check that we got ourselves a real number type
    try:
        if abs(RR(0.5) - 0.5) > 0.001:
            RR = RealField(53)  # hardcode something
    except TypeError:
        RR = RealField(53)  # hardcode something
    return RR(sqrt(2 * pi)) * stddev


@dataclass
class NoiseDistribution:
    """
    All noise distributions are instances of this class.
    """

    stddev: float
    mean: float = 0
    hamming_fraction: float = 1.0  # Hamming weight divided by dimension.

    def __lt__(self, other):
        try:
            return self.stddev < other.stddev
        except AttributeError:
            return self.stddev < other

    def __le__(self, other):
        try:
            return self.stddev <= other.stddev
        except AttributeError:
            return self.stddev <= other

    def __repr__(self):
        return f"D(σ={self.stddev:.2f}, μ={self.mean:.2f})"

    def __hash__(self):
        return hash((self.stddev, self.mean))

    def __len__(self):
        if hasattr(self, "n"):
            return self.n
        else:
            return 0

    @staticmethod
    def DiscreteGaussian(stddev, mean=0):
        """
        A discrete Gaussian distribution with standard deviation `stddev`.
        """
        return NoiseDistribution(stddev=RR(stddev), mean=RR(mean))

    @staticmethod
    def DiscreteGaussianAlpha(alpha, q, mean=0):
        """
        A discrete Gaussian distribution with standard deviation αq/√(2π).
        """
        stddev = stddevf(alpha * q)
        return NoiseDistribution.DiscreteGaussian(stddev=RR(stddev), mean=RR(mean))

    @staticmethod
    def CentredBinomial(eta):
        """
        Sample a_1, …, a_η, b_1, …, b_η and return Σ(a_i - b_i).
        """
        stddev = sqrt(eta / 2.0)
        # TODO: hamming_fraction
        return NoiseDistribution(stddev=RR(stddev), mean=RR(0))

    @staticmethod
    def Uniform(a, b):
        "Uniform distribution ∈ [a,b], endpoints inclusive."
        if b < a:
            raise ValueError(f"upper limit must be larger than lower limit but got: {b} < {a}")
        n = b - a + 1
        mean = (a + b) / RR(2)
        stddev = sqrt((n ** 2 - 1) / RR(12))

        if a <= 0 and 0 <= b:
            hamming_fraction = 1.0 / n
        else:
            hamming_fraction = 0.0

        return NoiseDistribution(stddev=stddev, mean=mean, hamming_fraction=hamming_fraction)

    @staticmethod
    def UniformMod(q):
        a = -(q // 2)
        b = q // 2
        if q % 2 == 0:
            b -= 1
        return NoiseDistribution.Uniform(a, b)

    @staticmethod
    def SparseTernary(n, p, m=None):
        if m is None:
            m = p
        mean = RR(p / n - m / n)
        stddev = RR(sqrt((p + m) / n))
        hamming_fraction = RR((p + m) / n)
        D = NoiseDistribution(stddev=stddev, mean=mean, hamming_fraction=hamming_fraction)
        D.h = p + m
        D.n = n
        return D
