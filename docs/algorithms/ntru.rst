.. _NTRU Primal Attacks:

NTRU Primal Attacks
=====================

We construct an (easy) example NTRU instance::

    from estimator import *
    params = NTRU.Parameters(n=200, q=7981, Xs=ND.UniformMod(3), Xe=ND.UniformMod(3))
    params

The simplest (and quickest to estimate) model is solving via uSVP and assuming the Geometric Series
Assumption (GSA) [Schnorr03]_. The success condition was formulated in [USENIX:ADPS16]_ and
studied/verified in [AC:AGVW17]_, [C:DDGR20]_, [PKC:PosVir21]_. The treatment of small secrets is
from [ACISP:BaiGal14]_. Here, the NTRU instance is treated as a homoegeneous LWE instance::

    NTRU.primal_usvp(params, red_shape_model="gsa")

We get a similar result if we use the ``GSA`` simulator. We do not get the identical result because
we optimize β and d separately::

    NTRU.primal_usvp(params, red_shape_model=Simulator.GSA)

To get a more precise answer we may use the CN11 simulator by Chen and Nguyen [AC:CheNgu11]_ (as
`implemented in FPyLLL
<https://github.com/fplll/fpylll/blob/master/src/fpylll/tools/bkz_simulator.py>_`)::

    NTRU.primal_usvp(params, red_shape_model=Simulator.CN11)

We can then improve on this result by first preprocessing the basis with block size β followed by a
single SVP call in dimension η [RSA:LiuNgu13]_. We call this the BDD approach since this is
essentially the same strategy as preprocessing a basis and then running a CVP solver::

    NTRU.primal_bdd(params, red_shape_model=Simulator.CN11)

We can improve these results further by exploiting the sparse secret in the hybrid attack
[C:HowgraveGraham07]_ guessing ζ positions of the secret::

    NTRU.primal_hybrid(params, red_shape_model=Simulator.CN11)

In addition to the primal secret key recovery attack, this module supports the dense sublattice
attack as formulated in [EC:KirFou17]_, and refined/verified in [AC:DucWoe21]_. The baseline
dense sublattice attack uses a 'z-shape' variant of the Geometric Series Assumption, called the
ZGSA::

    params.possibly_overstretched
    NTRU.primal_dsd(params, red_shape_model=Simulator.ZGSA)

Of course we can also use the CN11 simulator for this attack as well::

    NTRU.primal_dsd(params, red_dhape_model=Simulator.CN11)

**Note:** Currently, dense sublattice attack estimation is only supported if the distributions of
``f`` and ``g`` are equal. ``NTRU.primal_dsd()`` will return a ``NotImplementedError`` if this is
not the case. 
                        
