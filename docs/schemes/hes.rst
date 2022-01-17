Homomorphic Encryption Standard
===============================

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv111024128error
    >>> HESv111024128error
    LWEParameters(n=1024, q=134217728, Xs=D(σ=3.00), Xe=D(σ=3.00), m=1024, tag='HESv11error')
    >>> LWE.primal_bdd(HESv111024128error)
    rop: ≈2^143.3, red: ≈2^143.0, svp: ≈2^140.9, β: 374, η: 405, d: 2047, tag: bdd

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv111024128ternary
    >>> HESv111024128ternary
    LWEParameters(n=1024, q=134217728, Xs=D(σ=0.82), Xe=D(σ=3.00), m=1024, tag='HESv11ternary')
    >>> LWE.primal_hybrid(HESv111024128ternary)
    rop: ≈2^184.1, red: ≈2^183.3, svp: ≈2^182.9, β: 345, η: 2, ζ: 142, |S|: ≈2^225.1, d: 1907, prob: ≈2^-46.4, ...
   
