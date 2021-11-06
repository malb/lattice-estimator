.. _LWE Primal Attacks:

LWE Primal Attacks
==================
.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Flwe-primal.ipynb

We construct an (easy) example LWE instance::

    >>> from estimator import *
    >>> params = LWEParameters(n=200, q=7981, Xs=ND.SparseTernary(384, 16), Xe=ND.CenteredBinomial(4))
    >>> params
    LWEParameters(n=200, q=7981, Xs=D(σ=0.29), Xe=D(σ=1.41), m=+Infinity, tag=None)

The simplest (and quickest to estimate) model is solving via uSVP and assuming the Geometric Series
Assumption (GSA)::

    >>> primal_usvp(params, red_shape_model="gsa")
    rop: ≈2^53.3, red: ≈2^53.3, δ: 1.010838, β: 70, d: 362, tag: usvp

We get a similar result if we use the ``GSA`` simulator. We do not get the identical result because
we optimize β and d separately::

    >>> primal_usvp(params, red_shape_model=Simulator.GSA)
    rop: ≈2^53.7, red: ≈2^53.7, δ: 1.010720, β: 72, d: 331, tag: usvp

To get a more precise answer we may use the CN11 simulator by Chen and Nguyen (as implemented in FPLLL)::

    >>> primal_usvp(params, red_shape_model=Simulator.CN11)
    rop: ≈2^53.6, red: ≈2^53.6, δ: 1.010779, β: 71, d: 358, tag: usvp

We can then improve on this result by first preprocessing the basis with blocksize β followed by a
single SVP call in dimension η. We call this the BDD approach since this is essentially the same
strategy as preprocessing a basis and then running a CVP solver::

    >>> primal_bdd(params, red_shape_model=Simulator.CN11)
    rop: ≈2^51.5, red: ≈2^49.9, svp: ≈2^51.0, β: 58, η: 91, d: 345, tag: bdd

We can improve these results further by exploiting the sparse secret in the hybrid attack, guessing ζ
positions of the secret::

    >>> primal_hybrid(params, red_shape_model=Simulator.CN11)
    rop: ≈2^45.9, red: ≈2^45.5, svp: ≈2^43.8, β: 40, η: 2, ζ: 62, |S|: ≈2^43.4, d: 313, prob: 0.988, ...

