# -*- coding: utf-8 -*-
from sage.all import binomial, ZZ, log, ceil, RealField, oo
from sage.all import RealDistribution, RR, sqrt, prod


def babai(r, norm):
    """
    Babai probability following [Wun16]_.
    """
    R = [RR(sqrt(t) / (2 * norm)) for t in r]
    T = RealDistribution("beta", ((len(r) - 1) / 2, 1.0 / 2))
    probs = [1 - T.cum_distribution_function(1 - s ** 2) for s in R]
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
    :param majority: if `True` amplify a deicsional problem, not a computational one
       if `False` then we assume that we can check solutions, so one success suffices

    :returns: number of required trials to amplify
    """
    if target_success_probability < success_probability:
        return ZZ(1)
    if success_probability == 0.0:
        return oo

    prec = max(
        53,
        2 * ceil(abs(log(success_probability, 2))),
        2 * ceil(abs(log(1 - success_probability, 2))),
        2 * ceil(abs(log(target_success_probability, 2))),
        2 * ceil(abs(log(1 - target_success_probability, 2))),
    )
    prec = min(prec, 2048)
    RR = RealField(prec)

    success_probability = RR(success_probability)
    target_success_probability = RR(target_success_probability)

    try:
        if majority:
            eps = success_probability / 2
            return ceil(2 * log(2 - 2 * target_success_probability) / log(1 - 4 * eps ** 2))
        else:
            # target_success_probability = 1 - (1-success_probability)^trials
            return ceil(log(1 - target_success_probability) / log(1 - success_probability))
    except ValueError:
        return oo
