.. _LWE Dual Attacks:

LWE Dual Attacks
==================

We construct an (easy) example LWE instance::

    from estimator import *
    from estimator.lwe_dual import dual_hybrid, matzov
    params = LWE.Parameters(n=200, q=7981, Xs=ND.SparseTernary(16), Xe=ND.CenteredBinomial(4))
    params

The simplest (and quickest to estimate) algorithm is the "plain" dual attack as described in [PQCBook:MicReg09]_::

    LWE.dual(params)

We can improve these results by considering a dual hybrid attack as in [EC:Albrecht17]_, [INDOCRYPT:EspJouKha20]_::

    dual_hybrid(params)

Further improvements are possible using a meet-in-the-middle approach [IEEE:CHHS19]_::

   dual_hybrid(params, mitm_optimization=True)

We consider the variant from [MATZOV22]_::

   matzov(params)
