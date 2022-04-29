Homomorphic Encryption Standard
===============================

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv111024128error
    >>> HESv111024128error
    LWEParameters(n=1024, q=134217728, Xs=D(σ=3.00), Xe=D(σ=3.00), m=1024, tag='HESv11error')
    >>> LWE.primal_bdd(HESv111024128error)
    rop: ≈2^136.7, red: ≈2^136.4, svp: ≈2^134.3, β: 374, η: 405, d: 2047, tag: bdd

::

    >>> from estimator import *
    >>> from estimator.schemes import HESv111024128ternary
    >>> HESv111024128ternary
    LWEParameters(n=1024, q=134217728, Xs=D(σ=0.82), Xe=D(σ=3.00), m=1024, tag='HESv11ternary')
    >>> LWE.primal_hybrid(HESv111024128ternary)
    rop: ≈2^182.5, red: ≈2^181.7, svp: ≈2^181.4, β: 345, η: 2, ζ: 134, |S|: ≈2^212.4, d: 1915, prob: ≈2^-51.2, ↻: ≈2^53.4, tag: hybrid
   
