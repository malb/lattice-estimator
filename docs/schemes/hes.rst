Homomorphic Encryption Standard
===============================

::

   >>> from estimator import *
   >>> from estimator.schemes import HESv111024128error
   >>> HESv111024128error
   LWEParameters(n=1024, q=134217728, Xs=D(σ=3.00), Xe=D(σ=3.00), m=1024, tag='HESv11error')
   >>> primal_bdd(HESv111024128error)
   rop: ≈2^140.8, red: ≈2^140.5, svp: ≈2^138.5, β: 374, η: 405, d: 2047, tag: bdd

::

   >>> from estimator import *
   >>> from estimator.schemes import HESv111024128ternary
   >>> HESv111024128ternary
   LWEParameters(n=1024, q=134217728, Xs=D(σ=0.82), Xe=D(σ=3.00), m=1024, tag='HESv11ternary')
   >>> primal_hybrid(HESv111024128ternary)
   rop: ≈2^181.1, red: ≈2^180.2, svp: ≈2^179.9, β: 345, η: 2, ζ: 139, |S|: ≈2^220.3, ...
   
