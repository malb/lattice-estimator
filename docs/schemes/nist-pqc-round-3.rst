NIST PQC Round 3 Finalists
==========================

`Kyber <https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf>`__

::

   >>> from estimator import *
   >>> Kyber512
   LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.00), m=1024, tag='Kyber 512')
   >>> primal_bdd(Kyber512)
   rop: ≈2^137.8, red: ≈2^136.5, svp: ≈2^137.1, β: 365, η: 400, d: 981, tag: bdd

::

   >>> from estimator import *
   >>> Kyber768
   LWEParameters(n=768, q=3329, Xs=D(σ=1.00), Xe=D(σ=1.00), m=768, tag='Kyber 768')
   >>> primal_bdd(Kyber768)
   rop: ≈2^205.8, red: ≈2^204.3, svp: ≈2^205.1, β: 605, η: 642, d: 1423, tag: bdd

::

   >>> from estimator import *
   >>> Kyber1024
   LWEParameters(n=1024, q=3329, Xs=D(σ=1.00), Xe=D(σ=1.00), m=1024, tag='Kyber 1024')
   >>> primal_bdd(Kyber1024)
   rop: ≈2^276.1, red: ≈2^274.9, svp: ≈2^275.3, β: 854, η: 891, d: 1869, tag: bdd

`Saber <https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf>`__

::

   >>> from estimator import *
   >>> LightSaber
   LWEParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=2.00), m=512, tag='LightSaber')
   >>> primal_bdd(LightSaber)
   rop: ≈2^142.0, red: ≈2^141.1, svp: ≈2^141.0, β: 381, η: 414, d: 1010, tag: bdd

::

   >>> from estimator import *
   >>> Saber
   LWEParameters(n=768, q=8192, Xs=D(σ=1.41), Xe=D(σ=2.00), m=768, tag='Saber')
   >>> primal_bdd(Saber)
   rop: ≈2^209.7, red: ≈2^208.3, svp: ≈2^209.0, β: 619, η: 656, d: 1475, tag: bdd

::

   >>> from estimator import *
   >>> FireSaber
   LWEParameters(n=1024, q=8192, Xs=D(σ=1.22), Xe=D(σ=2.00), m=1024, tag='FireSaber')
   >>> primal_bdd(FireSaber)
   rop: ≈2^277.1, red: ≈2^275.7, svp: ≈2^276.5, β: 857, η: 895, d: 1884, tag: bdd


`NTRU <https://ntru.org/f/ntru-20190330.pdf>`__

::

   >>> from estimator import *
   >>> NTRUHPS2048509Enc
   LWEParameters(n=508, q=2048, Xs=D(σ=0.82), Xe=D(σ=0.71), m=508, tag='NTRUHPS2048509Enc')
   >>> primal_bdd(NTRUHPS2048509Enc)
   rop: ≈2^135.3, red: ≈2^133.9, svp: ≈2^134.6, β: 356, η: 391, d: 916, tag: bdd

::

   >>> from estimator import *
   >>> NTRUHPS2048677Enc
   LWEParameters(n=676, q=2048, Xs=D(σ=0.82), Xe=D(σ=0.61), m=676, tag='NTRUHPS2048677Enc')
   >>> primal_bdd(NTRUHPS2048677Enc)
   rop: ≈2^175.2, red: ≈2^174.0, svp: ≈2^174.4, β: 498, η: 533, d: 1179, tag: bdd

::

   >>> from estimator import *
   >>> NTRUHPS4096821Enc
   LWEParameters(n=820, q=4096, Xs=D(σ=0.82), Xe=D(σ=0.79), m=820, tag='NTRUHPS4096821Enc')
   >>> primal_bdd(NTRUHPS4096821Enc)
   rop: ≈2^204.4, red: ≈2^203.0, svp: ≈2^203.7, β: 600, η: 637, d: 1486, tag: bdd

::

   >>> from estimator import *
   >>> NTRUHRSS701Enc
   LWEParameters(n=700, q=8192, Xs=D(σ=0.82), Xe=D(σ=0.82), m=700, tag='NTRUHRSS701')
   >>> primal_bdd(NTRUHRSS701Enc)
   rop: ≈2^163.3, red: ≈2^162.5, svp: ≈2^162.0, β: 456, η: 489, d: 1288, tag: bdd
