.. _LWE Primal Attacks:

LWE Primal Attacks
==================

We construct an (easy) example LWE instance::

    from estimator import *
    params = LWE.Parameters(n=200, q=7981, Xs=ND.SparseTernary(384, 16), Xe=ND.CenteredBinomial(4))
    params

The simplest (and quickest to estimate) model is solving via uSVP and assuming the Geometric Series
Assumption (GSA) [Schnorr03]_. The success condition was formulated in [USENIX:ADPS16]_ and
studied/verified in [AC:AGVW17]_, [C:DDGR20]_, [PKC:PosVir21]_. The treatment of small secrets is
from [ACISP:BaiGal14]_::

    LWE.primal_usvp(params, red_shape_model="gsa")

We get a similar result if we use the ``GSA`` simulator. We do not get the identical result because
we optimize β and d separately::

    LWE.primal_usvp(params, red_shape_model=Simulator.GSA)

To get a more precise answer we may use the CN11 simulator by Chen and Nguyen [AC:CheNgu11]_ (as
`implemented in FPyLLL
<https://github.com/fplll/fpylll/blob/master/src/fpylll/tools/bkz_simulator.py>_`)::

    LWE.primal_usvp(params, red_shape_model=Simulator.CN11)

We can then improve on this result by first preprocessing the basis with block size β followed by a
single SVP call in dimension η [RSA:LiuNgu13]_. We call this the BDD approach since this is
essentially the same strategy as preprocessing a basis and then running a CVP solver::

    LWE.primal_bdd(params, red_shape_model=Simulator.CN11)

We can improve these results further by exploiting the sparse secret in the hybrid attack
[C:HowgraveGraham07]_ guessing ζ positions of the secret::

    LWE.primal_hybrid(params, red_shape_model=Simulator.CN11)
                        
