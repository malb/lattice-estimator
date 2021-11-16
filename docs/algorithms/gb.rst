.. _Arora-GB:

Arora-GB
========

We construct an (easy) example LWE instance and estimate the cost of solving it using Gröbner bases as described in [ICALP:AroGe11]_, [EPRINT:ACFP14]_::

    from estimator import *
    params = LWEParameters(n=64, q=7681, Xs=ND.DiscreteGaussian(3.0), Xe=ND.DiscreteGaussian(3.0), m=2^50)
    arora_gb(params)
    
The cost of this approach – Arora-GB – depends on the number of samples::

    arora_gb(params.updated(m=2^120))

If the noise distribution is bounded, this bounds the absolute degree and thus cost::

    arora_gb(params.updated(Xe=ND.UniformMod(7)))

Centered binomial distributions are also bounded::

    arora_gb(params.updated(Xe=ND.CenteredBinomial(8)))

The secret plays its role, too, in reducing the cost of solving::

    arora_gb(params.updated(Xs=ND.UniformMod(5), Xe=ND.CenteredBinomial(4), m=1024))

::

   arora_gb(params.updated(Xs=ND.UniformMod(3), Xe=ND.CenteredBinomial(4), m=1024))

