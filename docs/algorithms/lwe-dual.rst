.. _LWE Primal Attacks:

LWE Dual Attacks
==================

We construct an (easy) example LWE instance::

    from estimator import *
    params = LWEParameters(n=200, q=7981, Xs=ND.SparseTernary(384, 16), Xe=ND.CenteredBinomial(4))
    params

