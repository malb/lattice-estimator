.. _SIS Small-q Attacks:

SIS Small-q Attacks
===================

We construct an example SIS instance with a small modulus::

    from estimator import *
    params = SIS.Parameters(n=512, q=257, length_bound=801, m=1024, norm=2)
    params

The small-q attack from [C:DEP23]_ exploits the Z-shape geometry of BKZ-reduced bases for q-ary lattices when the modulus q is small relative to the length bound. The attack sieves in a projected sublattice and lifts short vectors over q-vectors::

    from estimator.sis_small_q import small_q
    small_q(params)

For on-the-fly lifting, which considers vectors during sieving rather than just the terminal database::

    small_q(params, otf_lift=True)

For the inhomogeneous ISIS problem with the optimized reduction from [C:DEP23]_::

    small_q(params, inhom="specific")

The attack is automatically included in ``SIS.estimate()`` when applicable::

    SIS.estimate(params)

