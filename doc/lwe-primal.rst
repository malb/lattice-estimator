.. _LWE Primal Attacks:

LWE Primal Attacks
==================

We construct an (easy) example LWE instance::

    from estimator import *
    params = LWEParameters(n=200, q=7981, Xs=ND.SparseTernary(384, 16), Xe=ND.CenteredBinomial(4))
    params

The simplest (and quickest to estimate) model is solving via uSVP and assuming the Geometric Series
Assumption (GSA)::

    primal_usvp(params, red_shape_model="gsa")

We get a similar result if we use the ``GSA`` simulator. We do not get the identical result because
we optimize β and d separately::

    primal_usvp(params, red_shape_model=Simulator.GSA)

To get a more precise answer we may use the CN11 simulator by Chen and Nguyen (as implemented in FPLLL)::

    primal_usvp(params, red_shape_model=Simulator.CN11)

We can then improve on this result by first preprocessing the basis with blocksize β followed by a
single SVP call in dimension η. We call this the BDD approach since this is essentially the same
strategy as preprocessing a basis and then running a CVP solver::

    primal_bdd(params, red_shape_model=Simulator.CN11)

We can improve these results further by exploiting the sparse secret in the hybrid attack, guessing ζ
positions of the secret::

    primal_hybrid(params, red_shape_model=Simulator.CN11)

