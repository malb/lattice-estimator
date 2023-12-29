# -*- coding: utf-8 -*-
from sage.all import binomial, ZZ, log, ceil, RealField, oo, exp, pi
from sage.all import RealDistribution, RR, sqrt, prod, erf
from .nd import sigmaf
from .conf import max_n_cache


chisquared_table = {i: None for i in range(2*max_n_cache+1)}
for i in range(2*max_n_cache+1):
    chisquared_table[i] = RealDistribution('chisquared', i)


def conditional_chi_squared(d1, d2, lt, l2):
    """
    Probability that a gaussian sample (var=1) of dim d1+d2 has length at most
    lt knowing that the d2 first coordinates have length at most l2

    :param d1: Dimension of non length-bounded coordinates
    :param d2: Dimension of length-bounded coordinates
    :param lt: Length threshold (maximum length of whole vector)
    :param l2: Length threshold for the first d2 coordinates.

    EXAMPLE::
        >>> from estimator import prob
        >>> prob.conditional_chi_squared(100, 5, 105, 1)
        0.6358492948586715

        >>> prob.conditional_chi_squared(100, 5, 105, 5)
        0.5764336909205551

        >>> prob.conditional_chi_squared(100, 5, 105, 10)
        0.5351747076352109

        >>> prob.conditional_chi_squared(100, 5, 50, 10)
        1.1707597206287592e-06

        >>> prob.conditional_chi_squared(100, 5, 50, .7)
        5.4021875103989546e-06
    """
    D1 = chisquared_table[d1].cum_distribution_function
    D2 = chisquared_table[d2].cum_distribution_function
    l2 = RR(l2)

    PE2 = D2(l2)
    # In large dim, we can get underflow leading to NaN
    # When this happens, assume lifting is successfully (underestimating security)
    if PE2==0:
        raise ValueError("Numerical underflow in conditional_chi_squared")

    steps = 5 * (d1 + d2)

    # Numerical computation of the integral
    proba = 0.
    for i in range(steps)[::-1]:
        l2_min = i * l2 / steps
        l2_mid = (i + .5) * l2 / steps
        l2_max = (i + 1) * l2 / steps

        PC2 = (D2(l2_max) - D2(l2_min)) / PE2
        PE1 = D1(lt - l2_mid)

        proba += PC2 * PE1

    return proba

def gaussian_cdf(mu, sigma, t):
    """
    Compute the cdf of a continuous gaussian random variable with mean mu and standard deviation
    sigma (i.e. computes Pr(X <= t), where X is a gaussian random variable).

    :params mu: the mean of the gaussian random variable.
    :params sigma: the standard deviation of the gaussian random variable.
    :params t: the limit at which to calculate the cdf.

    :returns: the evaluation of the cdf at t.
    """
    return RR((1/2)*(1 + erf((t - mu)/(sqrt(2)*sigma)))) 

def build_Gaussian_law(sigma, t):
    """
    Compute the pdf of a discrete gaussian of parameter sigma over the range [-t, t]

    :params sigma: the standard deviation of the gaussian random variable
    :params t: the (symmetric) limit of the range of the discrete gaussian.

    :returns: a dictionary of the pdf of the specified discrete gaussian, where D[i] = pr(X == i)
    """
    D = {}
    for i in range(0, t + 1):
        D[i] = exp(-i ** 2 / (2 * sigma ** 2))
        D[-i] = D[i]

    normalization = sum([D[i] for i in D])
    for i in D:
        D[i] = D[i] / normalization

    assert abs(sum([D[i] for i in range(-t, t + 1)]) - 1.) <= 10 ** -10
    return D


def mitm_babai_probability(r, stddev, q, fast=False):
    """
    Compute the "e-admissibility" probability associated to the mitm step, according to
    [EPRINT:SonChe19]_

    :params r: the squared GSO lengths
    :params stddev: the std.dev of the error distribution
    :params q: the LWE modulus
    :param fast: toggle for setting p = 1 (faster, but underestimates security)
    :return: probability for the mitm process

    # NOTE: the model sometimes outputs negative probabilities, we set p = 0 in this case
    """

    if fast:
        # overestimate the probability -> underestimate security
        return 1

    # get non-squared norms
    alphaq = sigmaf(stddev)
    probs = (
        RR(
            erf(s * sqrt(RR(pi)) / alphaq)
            + (alphaq / s) * ((exp(-s * sqrt(RR(pi)) / alphaq) - 1) / RR(pi))
        )
        for s in map(sqrt, r)
    )
    p = RR(prod(probs))
    return p if 0 <= p <= 1 else 0.0


def babai(r, norm):
    """
    Babai probability following [EPRINT:Wun16]_.

    """
    denom = float(2 * norm) ** 2
    T = RealDistribution("beta", ((len(r) - 1) / 2, 1.0 / 2))
    probs = [1 - T.cum_distribution_function(1 - r_ / denom) for r_ in r]
    return prod(probs)


def drop(n, h, k, fail=0, rotations=False):
    """
    Probability that ``k`` randomly sampled components have ``fail`` non-zero components amongst
    them.

    :param n: LWE dimension `n > 0`
    :param h: number of non-zero components
    :param k: number of components to ignore
    :param fail: we tolerate ``fail`` number of non-zero components amongst the `k` ignored
        components
    :param rotations: consider rotations of the basis to exploit ring structure (NTRU only)
    """

    N = n  # population size
    K = n - h  # number of success states in the population
    n = k  # number of draws
    k = n - fail  # number of observed successes
    prob_drop = binomial(K, k) * binomial(N - K, n - k) / binomial(N, n)
    if rotations:
        return 1 - (1 - prob_drop) ** N
    else:
        return prob_drop


def amplify(target_success_probability, success_probability, majority=False):
    """
    Return the number of trials needed to amplify current `success_probability` to
    `target_success_probability`

    :param target_success_probability: targeted success probability < 1
    :param success_probability: targeted success probability < 1
    :param majority: if `True` amplify a decisional problem, not a computational one
       if `False` then we assume that we can check solutions, so one success suffices

    :returns: number of required trials to amplify
    """
    if target_success_probability < success_probability:
        return ZZ(1)
    if success_probability == 0.0:
        return oo

    prec = max(
        53,
        2 * ceil(abs(float(log(success_probability, 2)))),
        2 * ceil(abs(float(log(1 - success_probability, 2)))),
        2 * ceil(abs(float(log(target_success_probability, 2)))),
        2 * ceil(abs(float(log(1 - target_success_probability, 2)))),
    )
    prec = min(prec, 2048)
    RR = RealField(prec)

    success_probability = RR(success_probability)
    target_success_probability = RR(target_success_probability)

    try:
        if majority:
            eps = success_probability / 2
            return ceil(2 * log(2 - 2 * target_success_probability) / log(1 - 4 * eps**2))
        else:
            # target_success_probability = 1 - (1-success_probability)^trials
            return ceil(log(1 - target_success_probability) / log(1 - success_probability))
    except ValueError:
        return oo


def amplify_sigma(target_advantage, sigma, q):
    """
    Amplify distinguishing advantage for a given σ and q

    :param target_advantage:
    :param sigma: (Lists of) Gaussian width parameters
    :param q: Modulus q > 0

    """
    try:
        sigma = sum(sigma_**2 for sigma_ in sigma).sqrt()
    except TypeError:
        pass
    if sigma > 16 * q:
        return oo

    advantage = float(exp(-float(pi) * (float(sigma / q) ** 2)))
    return amplify(target_advantage, advantage, majority=True)
