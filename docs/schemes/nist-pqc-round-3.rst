NIST PQC Round 3 Finalists
==========================

`Kyber <https://pq-crystals.org/kyber/data/kyber-specification-round3-20210804.pdf>`__

::

    >>> from estimator import *
    >>> schemes.Kyber512
    LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.22), m=512, tag='Kyber 512')
    >>> LWE.primal_bdd(schemes.Kyber512)
    rop: ≈2^140.3, red: ≈2^139.7, svp: ≈2^138.7, β: 391, η: 421, d: 1013, tag: bdd

::

    >>> from estimator import *
    >>> schemes.Kyber768
    LWEParameters(n=768, q=3329, Xs=D(σ=1.00), Xe=D(σ=1.00), m=768, tag='Kyber 768')
    >>> LWE.primal_bdd(schemes.Kyber768)
    rop: ≈2^201.0, red: ≈2^199.9, svp: ≈2^200.0, β: 606, η: 641, d: 1425, tag: bdd

::

    >>> from estimator import *
    >>> schemes.Kyber1024
    LWEParameters(n=1024, q=3329, Xs=D(σ=1.00), Xe=D(σ=1.00), m=1024, tag='Kyber 1024')
    >>> LWE.primal_bdd(schemes.Kyber1024)
    rop: ≈2^270.7, red: ≈2^269.8, svp: ≈2^269.6, β: 855, η: 890, d: 1873, tag: bdd

`Saber <https://www.esat.kuleuven.be/cosic/pqcrypto/saber/files/saberspecround3.pdf>`__

::

    >>> from estimator import *
    >>> schemes.LightSaber
    LWEParameters(n=512, q=8192, Xs=D(σ=1.58), Xe=D(σ=2.29, μ=-0.50), m=512, tag='LightSaber')
    >>> LWE.primal_bdd(schemes.LightSaber)
    rop: ≈2^140.0, red: ≈2^139.4, svp: ≈2^138.5, β: 390, η: 420, d: 1025, tag: bdd

::

    >>> from estimator import *
    >>> schemes.Saber
    LWEParameters(n=768, q=8192, Xs=D(σ=1.41), Xe=D(σ=2.29, μ=-0.50), m=768, tag='Saber')
    >>> LWE.primal_bdd(schemes.Saber)
    rop: ≈2^208.1, red: ≈2^206.9, svp: ≈2^207.3, β: 631, η: 667, d: 1478, tag: bdd

::

    >>> from estimator import *
    >>> schemes.FireSaber
    LWEParameters(n=1024, q=8192, Xs=D(σ=1.22), Xe=D(σ=2.29, μ=-0.50), m=1024, tag='FireSaber')
    >>> LWE.primal_bdd(schemes.FireSaber)
    rop: ≈2^275.8, red: ≈2^274.9, svp: ≈2^274.7, β: 873, η: 908, d: 1894, tag: bdd


`NTRU <https://ntru.org/f/ntru-20190330.pdf>`__

::

    >>> from estimator import *
    >>> schemes.NTRUHPS2048509Enc
    NTRUParameters(n=508, q=2048, Xs=D(σ=0.82), Xe=T(p=127, m=127, n=508), m=508, tag='NTRUHPS2048509Enc', ntru_type='matrix')
    >>> NTRU.primal_bdd(schemes.NTRUHPS2048509Enc)
    rop: ≈2^131.1, red: ≈2^130.1, svp: ≈2^130.1, β: 357, η: 390, d: 916, tag: bdd

::

    >>> from estimator import *
    >>> schemes.NTRUHPS2048677Enc
    NTRUParameters(n=676, q=2048, Xs=D(σ=0.82), Xe=T(p=127, m=127, n=676), m=676, tag='NTRUHPS2048677Enc', ntru_type='matrix')
    >>> NTRU.primal_bdd(schemes.NTRUHPS2048677Enc)
    rop: ≈2^170.7, red: ≈2^169.6, svp: ≈2^169.9, β: 498, η: 533, d: 1179, tag: bdd

::

    >>> from estimator import *
    >>> schemes.NTRUHPS4096821Enc
    NTRUParameters(n=820, q=4096, Xs=D(σ=0.82), Xe=T(p=255, m=255, n=820), m=820, tag='NTRUHPS4096821Enc', ntru_type='matrix')
    >>> NTRU.primal_bdd(schemes.NTRUHPS4096821Enc)
    rop: ≈2^199.6, red: ≈2^198.6, svp: ≈2^198.6, β: 601, η: 636, d: 1485, tag: bdd

::

    >>> from estimator import *
    >>> schemes.NTRUHRSS701Enc
    NTRUParameters(n=700, q=8192, Xs=D(σ=0.82), Xe=D(σ=0.82), m=700, tag='NTRUHRSS701', ntru_type='matrix')
    >>> NTRU.primal_bdd(schemes.NTRUHRSS701Enc)
    rop: ≈2^158.6, red: ≈2^157.6, svp: ≈2^157.6, β: 454, η: 489, d: 1306, tag: bdd
