Homomorphic Encryption Parameters
=================================

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv111024128error
    >>> HESv111024128error
    LWEParameters(n=1024, q=134217728, Xs=D(σ=3.00), Xe=D(σ=3.00), m=1024, tag='HESv11error')
    >>> LWE.primal_bdd(HESv111024128error)
    rop: ≈2^136.7, red: ≈2^136.4, svp: ≈2^134.3, β: 374, η: 404, d: 2044, tag: bdd

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv111024128ternary
    >>> HESv111024128ternary
    LWEParameters(n=1024, q=134217728, Xs=D(σ=0.82), Xe=D(σ=3.00), m=1024, tag='HESv11ternary')
    >>> LWE.primal_hybrid(HESv111024128ternary)
    rop: ≈2^183.5, red: ≈2^182.6, svp: ≈2^182.4, β: 345, η: 2, ζ: 134, |S|: ≈2^212.4, d: 1914, prob: ≈2^-52.2, ↻: ≈2^54.4, tag: hybrid
   
::

    >>> from estimator import *
    >>> from estimator.schemes import SEAL22_8192
    >>> SEAL22_8192
    LWEParameters(n=8192, q=107839786668602559178668060348078522694548577690162289924414373888001, Xs=D(σ=0.82), Xe=D(σ=3.19), m=+Infinity, tag='SEAL22_8192')
    >>> LWE.dual_hybrid(SEAL22_8192)
    rop: ≈2^121.8, red: ≈2^121.8, guess: ≈2^101.7, β: 306, p: 3, ζ: 10, t: 40, β': 331, N: ≈2^68.1, m: ≈2^13.0
