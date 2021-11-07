Homomorphic Encryption Standard
===============================

::

   >>> from estimator import *
   >>> from estimator.schemes import HESv111024128error
   >>> HESv111024128error
   LWEParameters(n=1024, q=134217728, Xs=D(σ=3.00), Xe=D(σ=3.00), m=1024, tag='HESv11error')
   >>> primal_bdd(HESv111024128error)
   rop: ≈2^140.8, red: ≈2^140.2, svp: ≈2^139.3, β: 373, η: 408, d: 11, tag: bdd

::

   >>> from estimator import *
   >>> from estimator.schemes import HESv111024128ternary
   >>> HESv111024128ternary
   LWEParameters(n=1024, q=134217728, Xs=D(σ=0.82), Xe=D(σ=3.00), m=1024, tag='HESv11ternary')
   >>> primal_hybrid(HESv111024128ternary)
   rop: ≈2^130.6, red: ≈2^130.1, svp: ≈2^129.0, β: 337, η: 371, d: 1966, tag: bdd
   
