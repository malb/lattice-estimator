Homomorphic Encryption Parameters
=================================
Import a selection of parameters used in applications of Homomorphic Encryption, and estimate the best performing attack against each.

::

    >>> from estimator import *
    >>> from estimator.schemes import TFHErs_LWE
    >>> LWE.dual_hybrid(TFHErs_LWE)
    rop: ≈2^134.9, red: ≈2^134.9, guess: ≈2^128.0, β: 366, p: 2, ζ: 0, t: 110, β': 375, N: ≈2^77.1, m: 918

::

    >>> from estimator import *
    >>> from estimator.schemes import TFHErs_RLWE
    >>> LWE.dual_hybrid(TFHErs_RLWE)
    rop: ≈2^134.8, red: ≈2^134.7, guess: ≈2^128.0, β: 360, p: 2, ζ: 0, t: 110, β': 375, N: ≈2^77.3, m: ≈2^11.0

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv2_8192_128_ternary
    >>> LWE.primal_bdd(HESv2_8192_128_ternary)
    rop: ≈2^128.3, red: ≈2^128.0, svp: ≈2^125.7, β: 332, η: 373, d: 16216, tag: bdd

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv2_32768_128_ternary
    >>> LWE.primal_bdd(HESv2_32768_128_ternary)
    rop: ≈2^128.0, red: ≈2^128.0, svp: ≈2^118.7, β: 325, η: 348, d: 63105, tag: bdd

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv2_131072_128_ternary
    >>> LWE.primal_bdd(HESv2_131072_128_ternary)
    rop: ≈2^127.9, red: ≈2^127.9, svp: ≈2^118.7, β: 317, η: 348, d: 260026, tag: bdd
