# -*- coding: utf-8 -*-

from copy import copy
from dataclasses import dataclass

from sage.all import binomial, ceil, exp, floor, log, oo, parent, pi, QQ, RealField, RR, sqrt


def stddevf(sigma):
    """
    Gaussian width parameter σ → standard deviation.

    :param sigma: Gaussian width parameter σ

    EXAMPLE::

        >>> from estimator import *
        >>> ND.stddevf(64.0)
        25.532...

        >>> ND.stddevf(64)
        25.532...

        >>> ND.stddevf(RealField(256)(64)).prec()
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
    Standard deviation → Gaussian width parameter σ.

    :param stddev: standard deviation

    EXAMPLE::

        >>> from estimator import *
        >>> n = 64.0
        >>> ND.sigmaf(ND.stddevf(n))
        64.000...

        >>> ND.sigmaf(RealField(128)(1.0))
        2.5066282746310005024157652848110452530
        >>> ND.sigmaf(1.0)
        2.506628274631...
        >>> ND.sigmaf(1)
        2.506628274631...
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
    It is recommended to pick one of the following available implementations below:
    - DiscreteGaussian
    - DiscreteGaussianAlpha
    - CenteredBinomial
    - Uniform
    - UniformMod
    - SparseTernary
    - SparseBinary
    - Binary
    - Ternary

    NOTE:
    Generally, to generate an LWE parameter you call one of the above for the secret and error,
    **without** specifying the dimension `n` and `m` for secret/error respectively!
    These are initialized, when constructing the LWEParameters object.
    """
    n: int = None  # dimension of noise
    mean: float = 0  # expectation value
    stddev: float = 0  # standard deviation (square root of variance)
    bounds: tuple = (-oo, oo)  # range in which each coefficient is sampled with high probability
    is_Gaussian_like: bool = False  # whether the distribution "decays like a gaussian"
    _density: float = 1.0  # proportion of nonzero coefficients in a sample

    def __lt__(self, other):
        """
        We compare distributions by comparing their standard deviation.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussian(2.0) < ND.CenteredBinomial(18)
            True
            >>> ND.DiscreteGaussian(3.0) < ND.CenteredBinomial(18)
            False
            >>> ND.DiscreteGaussian(4.0) < ND.CenteredBinomial(18)
            False

        """
        try:
            return self.stddev < other.stddev
        except AttributeError:
            return self.stddev < other

    def __le__(self, other):
        """
        We compare distributions by comparing their standard deviation.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussian(2.0) <= ND.CenteredBinomial(18)
            True
            >>> ND.DiscreteGaussian(3.0) <= ND.CenteredBinomial(18)
            True
            >>> ND.DiscreteGaussian(4.0) <= ND.CenteredBinomial(18)
            False

        """
        try:
            return self.stddev <= other.stddev
        except AttributeError:
            return self.stddev <= other

    def __str__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussianAlpha(0.01, 7681)
            D(σ=30.64)

        """
        if self.n:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f}, n={int(self.n)})"
        else:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f})"

    def __repr__(self):
        if self.mean == 0.0:
            return f"D(σ={float(self.stddev):.2f})"
        else:
            return f"D(σ={float(self.stddev):.2f}, μ={float(self.mean):.2f})"

    def __hash__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> hash(ND.DiscreteGaussian(3.0, 1.0)) == hash((3.0, 1.0, None))
            True

        """
        return hash((self.stddev, self.mean, self.n))

    def __len__(self):
        """
        Dimension of this noise distribution, i.e. number of coefficients that gets sampled.

        EXAMPLE::

            >>> from estimator import *
            >>> len(ND.SparseTernary(128, n=1024))
            1024

        """
        if self.n is None:
            raise ValueError("Distribution has no length.")
        return self.n

    def resize(self, new_n):
        """
        Return an altered distribution having a dimension `new_n`.

        :param int new_n: new dimension to change to
        """
        new_self = copy(self)
        new_self.n = new_n
        return new_self

    @property
    def hamming_weight(self):
        """
        The number of non-zero coefficients in this distribution

        EXAMPLE::

            >>> from estimator import *
            >>> ND.SparseTernary(128, n=1024).hamming_weight
            256
            >>> ND.SparseTernary(128, 64, 1024).hamming_weight
            192
        """
        return round(len(self) * float(self._density))

    @property
    def is_bounded(self):
        """
        Whether the value of coefficients are bounded
        """
        return (self.bounds[1] - self.bounds[0]) < oo

    @property
    def is_sparse(self):
        """
        Whether the density of the distribution is < 1/2.
        Note: 1/2 might be considered somewhat arbitrary.
        """
        # NOTE: somewhat arbitrary
        return self._density < 0.5

    def support_size(self, fraction=1.0):
        raise NotImplementedError("support_size")


class DiscreteGaussian(NoiseDistribution):
    """
    A discrete Gaussian distribution with standard deviation ``stddev`` per component.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.DiscreteGaussian(3.0, 1.0)
        D(σ=3.00, μ=1.00)
    """
    # cut-off for Gaussian distributions
    gaussian_tail_bound: int = 2
    # probability that a coefficient falls within the cut-off
    gaussian_tail_prob: float = 1 - 2 * exp(-4 * pi)

    def __init__(self, stddev, mean=0, n=None):
        stddev, mean = RR(stddev), RR(mean)
        b_val = oo if n is None else ceil(log(n, 2) * stddev)
        density = max(0.0, 1 - RR(1 / sigmaf(stddev)))  # NOTE: approximation that is accurate for large stddev.

        super().__init__(
            n=n,
            mean=mean,
            stddev=stddev,
            bounds=(-b_val, b_val),
            _density=density,
            is_Gaussian_like=True,
        )

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.DiscreteGaussian(1.0, n=128).support_size(0.99)
            2.68643790357272e174
        """
        # We will treat this noise distribution as bounded with failure probability `1 - fraction`.
        n = len(self)
        t = self.gaussian_tail_bound
        p = self.gaussian_tail_prob

        if p**n < fraction:
            raise NotImplementedError(
                f"TODO(DiscreteGaussian.support_size): raise t. {RR(p ** n)}, {n}, {fraction}"
            )

        b = 2 * t * sigmaf(self.stddev) + 1
        return RR(2.0 * b + 1)**n


def DiscreteGaussianAlpha(alpha, q, mean=0, n=None):
    """
    A discrete Gaussian distribution with standard deviation α⋅q/√(2π) per component.

    EXAMPLE::

        >>> from estimator import *
        >>> alpha, q = 0.001, 2048
        >>> ND.DiscreteGaussianAlpha(alpha, q)
        D(σ=0.82)
        >>> ND.DiscreteGaussianAlpha(alpha, q) == ND.DiscreteGaussian(ND.stddevf(alpha * q))
        True
    """
    return DiscreteGaussian(stddevf(alpha * q), mean, n)


class CenteredBinomial(NoiseDistribution):
    """
    Sample a_1, …, a_η, b_1, …, b_η uniformly from {0, 1}, and return Σ(a_i - b_i).

    EXAMPLE::

        >>> from estimator import *
        >>> ND.CenteredBinomial(8)
        D(σ=2.00)
    """
    def __init__(self, eta, n=None):
        density = 1 - binomial(2 * eta, eta) * 2 ** (-2 * eta)

        super().__init__(
            n=n,
            mean=0,
            stddev=RR(sqrt(eta / 2.0)),
            bounds=(-eta, eta),
            _density=density,
            is_Gaussian_like=True,
        )

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.CenteredBinomial(3, 10).support_size()
            282475249
            >>> ND.CenteredBinomial(3, 10).support_size(0.99)
            279650497
        """
        # TODO: this might be suboptimal/inaccurate for binomial distribution
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**len(self))


class Uniform(NoiseDistribution):
    """
    Uniform distribution ∈ ``ZZ ∩ [a, b]``, endpoints inclusive.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.Uniform(-3, 3)
        D(σ=2.00)
        >>> ND.Uniform(-4, 3)
        D(σ=2.29, μ=-0.50)
    """
    def __init__(self, a, b, n=None):
        a, b = int(ceil(a)), int(floor(b))
        if b < a:
            raise ValueError(f"upper limit must be larger than lower limit but got: {b} < {a}")
        m = b - a + 1

        super().__init__(
            n=n,
            mean=RR((a + b) / 2),
            stddev=RR(sqrt((m**2 - 1) / 12)),
            bounds=(a, b),
            _density=(1 - 1 / m if a <= 0 and b >= 0 else 1),
        )

    def __hash__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> hash(ND.Uniform(-10, 10)) == hash(("Uniform", (-10, 10), None))
            True
        """
        return hash(("Uniform", self.bounds, self.n))

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.Uniform(-3, 3, 64).support_size(0.99)
            1207562882759477428726191443614714994252339953407098880
        """
        # TODO: this might be suboptimal/inaccurate for binomial distribution
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**len(self))


def UniformMod(q, n=None):
    """
    Uniform mod ``q``, with balanced representation, i.e. values in ZZ ∩ [-q/2, q/2).

    EXAMPLE::

        >>> from estimator import *
        >>> ND.UniformMod(7)
        D(σ=2.00)
        >>> ND.UniformMod(8)
        D(σ=2.29, μ=-0.50)
        >>> ND.UniformMod(2) == ND.Uniform(-1, 0)
        True
    """
    a = -(q // 2)
    b = a + q - 1
    return Uniform(a, b, n=n)


class TUniform(NoiseDistribution):
    """
    TUniform distribution ∈ ``ZZ ∩ [-2**b, 2**b]``, endpoints inclusive.
    This distribution samples the two end-points with probability 1/2**(b+2) and the
    intermediate points with probability 1/2**(b+1).

    As an example, with b=0 this distribution samples ±1 each with probability 1/4 and
    0 with probability 1/2.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.TUniform(0)
        D(σ=0.71)
        >>> ND.TUniform(10)
        D(σ=591.21)
    """
    def __init__(self, b, n=None):
        b = int(ceil(b))

        super().__init__(
            n=n,
            mean=RR(0),
            stddev=RR(sqrt((2**(2*b+1) + 1)/6)),
            bounds=(-2**b, 2**b),
            _density=(1 - 1 / 2**(b+1)),
        )

    def __hash__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> hash(ND.TUniform(2)) == hash(("TUniform", (-4, 4), None))
            True
        """
        return hash(("TUniform", self.bounds, self.n))

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.TUniform(0, 64).support_size(0.99)
            3399346982089587232319333203968
        """
        a, b = self.bounds
        return ceil(RR(fraction) * (b - a + 1)**len(self))


class SparseTernary(NoiseDistribution):
    """
    Distribution of vectors of length ``n`` with ``p`` entries of 1 and ``m`` entries of -1, rest 0.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.SparseTernary(10, n=100)
        T(p=10, m=10, n=100)
        >>> ND.SparseTernary(10, 10, 100)
        T(p=10, m=10, n=100)
        >>> ND.SparseTernary(10, 8, 100)
        T(p=10, m=8, n=100)
        >>> ND.SparseTernary(0, 0, 0).support_size()
        1
    """
    def __init__(self, p, m=None, n=None):
        p, m = int(p), int(p if m is None else m)
        self.p, self.m = p, m

        # Yes, n=0 might happen when estimating the cost of the dual attack! Support size is 1
        if n is None:
            # Treat it the same as n=0.
            n = 0
        mean = 0 if n == 0 else RR((p - m) / n)
        density = 0 if n == 0 else RR((p + m) / n)
        stddev = sqrt(density - mean**2)

        super().__init__(
            n=n,
            mean=mean,
            stddev=stddev,
            bounds=(0 if m == 0 else -1, 0 if p == 0 else 1),
            _density=density,
        )

    def __hash__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> hash(ND.SparseTernary(16, n=128)) == hash(("SparseTernary", 128, 16, 16))
            True
        """
        return hash(("SparseTernary", self.n, self.p, self.m))

    def resize(self, new_n):
        """
        Return an altered distribution having a dimension `new_n`.
        Assumes `p` and `m` stay the same.
        """
        return SparseTernary(self.p, self.m, new_n)

    def split_balanced(self, new_n, new_hw=None):
        """
        Split the +1 and -1 entries in a balanced way, and return 2 SparseTernary distributions:
        one of dimension `new_n` and the other of dimension `n - new_n`.

        :param new_n: dimension of the first noise distribution
        :param new_hw: hamming weight of the first noise distribution. If none, we take the most likely weight.
        :return: tuple of (SparseTernary, SparseTernary)
        """
        n, hw = len(self), self.hamming_weight
        if new_hw is None:
            # Most likely split has same density: new_hw / new_n = hw / n.
            new_hw = int(QQ(hw * new_n / n).round('down'))

        new_p = int((QQ(new_hw * self.p) / hw).round('down'))
        new_m = new_hw - new_p
        return (
            SparseTernary(new_p, new_m, new_n),
            SparseTernary(self.p - new_p, self.m - new_m, n - new_n)
        )

    def split_probability(self, new_n, new_hw=None):
        """
        Compute probability of splitting in a way that one half having `new_n` coefficients has
        `new_hw` of the weight, and the remaining part the rest. This is naturally the proportion
        of such splits divided this support size.
        """
        left, right = self.split_balanced(new_n, new_hw)
        return left.support_size() * right.support_size() / self.support_size()

    @property
    def is_sparse(self):
        """
        Always say this is a sparse distribution, even if p + m >= n/2, because there is correlation between the
        coefficients: if you split the distribution into two of half the length, then you expect in each of them to be
        half the weight.
        """
        return True

    @property
    def hamming_weight(self):
        return self.p + self.m

    def support_size(self, fraction=1.0):
        """
        Compute the size of the support covering the probability given as fraction.

        EXAMPLE::

            >>> from estimator import *
            >>> ND.SparseTernary(8, 8, 64).support_size()
            6287341680214194176
        """
        n, p, m = len(self), self.p, self.m
        return ceil(binomial(n, p) * binomial(n - p, m) * RR(fraction))

    def __str__(self):
        """
        EXAMPLE::

            >>> from estimator import *
            >>> ND.SparseTernary(20, 20, n=100)
            T(p=20, m=20, n=100)

        """
        if self.n:
            return f"T(p={self.p}, m={self.m}, n={int(self.n)})"
        else:
            return f"T(p={int(self.p)}, m={int(self.m)})"

    def __repr__(self):
        return str(self)


def SparseBinary(hw, n=None):
    """
    Sparse binary noise distribution having `hw` coefficients equal to 1, and the rest zero.

    EXAMPLE::

        >>> from estimator import *
        >>> ND.SparseBinary(10).bounds
        (0, 1)
    """
    return SparseTernary(hw, 0, n)


"""
Binary noise uniform from {0, 1}^n
"""
Binary = Uniform(0, 1)

"""
Ternary noise uniform from {-1, 0, 1}^n
"""
Ternary = Uniform(-1, 1)
