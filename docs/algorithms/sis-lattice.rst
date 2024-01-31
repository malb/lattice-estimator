.. _SIS Lattice Attacks:

SIS Lattice Attacks
=====================

We construct an (easy) example SIS instance::

    from estimator import *
    params = SIS.Parameters(n=113, q=2048, length_bound=512, norm="l2") 
    params

The simplest (and quickest to estimate) model is solving for the SIS instance with a euclidian norm length bound and assuming the Gaussian heuristic [CheNgu12]_.Then, we can solve for the required root hermite factor [EC:GamNgu08]_ that will guarantee BKZ outputs a short enough vector::

    SIS.lattice(params)

The exact reduction shape model doesn't matter when using euclidian norm bounds, as the required block size is calculated directly from the length bound. 

For infinity norm length bounds, we have two separate analyses. Both follow the same basic strategy. We use the worst case euclidian norm bound as a lower bound on the hardness. Then, we analyze the probability of obtaining a short vector where every coordinate meets the infinity norm constraint. When sqrt(m)*length_bound is less than the modulus q, we follow the analysis of the MATZOV report ([MATZOV22]_ P.18). We simulate the cost of generating *many* short vectors and treat each coordinate of the vector as an i.i.d Gaussian random variable with standard deviation equal to the length(s) of these short vectors divided by the square root of the dimension.::

    params = SIS.Parameters(n=113, q=2048, length_bound=50, norm="linf")
    SIS.lattice(params)

When sqrt(m)*length_bound is **greater than** the modulus, we follow the analysis present in the NIST round 3 Dilithium specification ([Dilithium21]_ P.35). Here, since BKZ can now produce q vectors at the given length bound (which will always satisfy the bound), we explicitly account for the q-ary structure of the lattice. Every coordinate corresponding to a q-vector yields uniformly random values, while the middle region of the basis produces Gaussian random variables as above. To explicitly account for this q-ary structure, use the ``ZGSA`` simulator.:: 

    SIS.lattice(params.updated(length_bound=70)), red_shape_model=Simulator.ZGSA)

To get a more precise answer we may use the CN11 simulator by Chen and Nguyen [AC:CheNgu11]_ (as `implemented in FPyLLL <https://github.com/fplll/fpylll/blob/master/src/fpylll/tools/bkz_simulator.py>`__)::

    SIS.lattice(params.updated(length_bound=70)), red_shape_model=Simulator.CN11)

Another option is to simulate a rerandomization of the basis, such that the q-vectors are *forgotten*. This results in the ``LGSA`` simulator, where the short, unit vectors are still present in the basis. See Figure 12 in the dilithium submission for an example.We can then improve on this result by first preprocessing the basis with block size β followed by a single SVP call in dimension η [RSA:LiuNgu13]_. We call this the BDD approach since this is essentially the same strategy as preprocessing a basis and then running a CVP solver::

    SIS.lattice(params.updated(length_bound=70)), red_shape_model=Simulator.LGSA)

**Note:** Currently, lattice attack estimation is only avalailable for euclidian (``l2``) and infinity (``linf``) norms. ``SIS.lattice()`` will return a ``NotImplementedError`` if one of these two norms are not selected.
                        
