NIST PQC Round 3 Finalists
==========================

`Kyber <https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf>`__

::

    >>> from estimator import *
    >>> Kyber512
    LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.22), m=512, tag='Kyber 512')
    >>> LWE.primal_bdd(Kyber512)
    rop: ≈2^146.9, red: ≈2^146.3, svp: ≈2^145.4, β: 391, η: 421, d: 1013, tag: bdd

::

    >>> from estimator import *
    >>> Kyber768
    LWEParameters(n=768, q=3329, Xs=D(σ=1.00), Xe=D(σ=1.00), m=768, tag='Kyber 768')
    >>> LWE.primal_bdd(Kyber768)
    rop: ≈2^208.2, red: ≈2^207.1, svp: ≈2^207.2, β: 606, η: 641, d: 1425, tag: bdd

::

    >>> from estimator import *
    >>> Kyber1024
    LWEParameters(n=1024, q=3329, Xs=D(σ=1.00), Xe=D(σ=1.00), m=1024, tag='Kyber 1024')
    >>> LWE.primal_bdd(Kyber1024)
    rop: ≈2^278.5, red: ≈2^277.6, svp: ≈2^277.5, β: 855, η: 890, d: 1873, tag: bdd

`Saber <https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf>`__

::

    >>> from estimator import *
    >>> LightSaber
    LWEParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=2.00), m=512, tag='LightSaber')
    >>> LWE.primal_bdd(LightSaber)
    rop: ≈2^144.4, red: ≈2^143.8, svp: ≈2^142.9, β: 382, η: 412, d: 1024, tag: bdd

::

    >>> from estimator import *
    >>> Saber
    LWEParameters(n=768, q=8192, Xs=D(σ=1.41), Xe=D(σ=2.00), m=768, tag='Saber')
    >>> LWE.primal_bdd(Saber)
    rop: ≈2^212.1, red: ≈2^211.1, svp: ≈2^211.2, β: 620, η: 655, d: 1475, tag: bdd

::

    >>> from estimator import *
    >>> FireSaber
    LWEParameters(n=1024, q=8192, Xs=D(σ=1.22), Xe=D(σ=2.00), m=1024, tag='FireSaber')
    >>> LWE.primal_bdd(FireSaber)
    rop: ≈2^279.5, red: ≈2^278.5, svp: ≈2^278.6, β: 858, η: 894, d: 1886, tag: bdd


`NTRU <https://ntru.org/f/ntru-20190330.pdf>`__

::

    >>> from estimator import *
    >>> NTRUHPS2048509Enc
    LWEParameters(n=508, q=2048, Xs=D(σ=0.82), Xe=D(σ=0.71), m=508, tag='NTRUHPS2048509Enc')
    >>> LWE.primal_bdd(NTRUHPS2048509Enc)
    rop: ≈2^137.7, red: ≈2^136.6, svp: ≈2^136.7, β: 357, η: 390, d: 916, tag: bdd

::

    >>> from estimator import *
    >>> NTRUHPS2048677Enc
    LWEParameters(n=676, q=2048, Xs=D(σ=0.82), Xe=D(σ=0.61), m=676, tag='NTRUHPS2048677Enc')
    >>> LWE.primal_bdd(NTRUHPS2048677Enc)
    rop: ≈2^177.6, red: ≈2^176.7, svp: ≈2^176.6, β: 499, η: 532, d: 1181, tag: bdd

::

    >>> from estimator import *
    >>> NTRUHPS4096821Enc
    LWEParameters(n=820, q=4096, Xs=D(σ=0.82), Xe=D(σ=0.79), m=820, tag='NTRUHPS4096821Enc')
    >>> LWE.primal_bdd(NTRUHPS4096821Enc)
    rop: ≈2^206.8, red: ≈2^205.8, svp: ≈2^205.8, β: 601, η: 636, d: 1485, tag: bdd

::

    >>> from estimator import *
    >>> NTRUHRSS701Enc
    LWEParameters(n=700, q=8192, Xs=D(σ=0.82), Xe=D(σ=0.82), m=700, tag='NTRUHRSS701')
    >>> LWE.primal_bdd(NTRUHRSS701Enc)
    rop: ≈2^165.7, red: ≈2^164.7, svp: ≈2^164.8, β: 455, η: 490, d: 1294, tag: bdd
