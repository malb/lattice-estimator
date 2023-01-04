.. _Coded-BKW for LWE:

Coded-BKW for LWE
=================

We construct an example LWE instance::

    from estimator import *
    params = LWE.Parameters(n=400, q=7981, Xs=ND.SparseTernary(384, 16), Xe=ND.CenteredBinomial(4), m=800)
    params

and estimate the cost of Coded-BKW [C:GuoJohSta15]_, [C:KirFou15]_::

    LWE.coded_bkw(params)

All BKW variants require *a lot* of samples. As called above, the algorithm will produce the required samples from what it is given, which increases the noise. Let's pretend we have as many as we like::

    LWE.coded_bkw(params.updated(m=2**100))
