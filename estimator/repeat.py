# -*- coding: utf-8 -*-
"""
Repeat algorithms.
"""

from sage.all import pi, exp, oo, ceil, log, RealField


def amplify(target_prob, prob, majority=False):
    """
    Return the number of trials needed to amplify current ``prob`` to ``target_prob``

    :param target_prob: Targeted success probability < 1.
    :param prob: Success probability < 1.
    :param majority: if ``True`` amplify a deicsional problem, not a computational one
       if ``False`` then we assume that we can check solutions, so one success suffices

    :returns: number of required trials to amplify.
    """
    if target_prob < prob:
        return int(1)
    if prob == 0.0:
        return oo

    prec = max(
        53,
        2 * ceil(abs(log(prob, 2))),
        2 * ceil(abs(log(1 - prob, 2))),
        2 * ceil(abs(log(target_prob, 2))),
        2 * ceil(abs(log(1 - target_prob, 2))),
    )
    prec = min(prec, 2048)
    RR = RealField(prec)

    prob = RR(prob)
    target_prob = RR(target_prob)

    try:
        if majority:
            eps = prob / 2
            return ceil(2 * log(2 - 2 * target_prob) / log(1 - 4 * eps ** 2))
        else:
            # target_prob = 1 - (1-prob)^trials
            return ceil(log(1 - target_prob) / log(1 - prob))
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
        sigma = sum(sigma_ ** 2 for sigma_ in sigma).sqrt()
    except TypeError:
        pass
    advantage = float(exp(-float(pi) * (float(sigma / q) ** 2)))
    return amplify(target_advantage, advantage, majority=True)
