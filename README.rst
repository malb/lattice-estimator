Security Estimates for Lattice Problems
=======================================

This `Sage <http://sagemath.org>`__ module provides functions for estimating the concrete security
of `Learning with Errors <https://en.wikipedia.org/wiki/Learning_with_errors>`__ instances.

The main intend of this estimator is to give designers an easy way to choose parameters resisting
known attacks and to enable cryptanalysts to compare their results and ideas with other techniques
known in the literature.

Usage Examples
--------------

::

    sage: from estimator import *
    sage: Kyber512
    LWEParameters(n=512, q=3329, Xs=D(σ=1.22, μ=0.00), Xe=D(σ=1.00, μ=0.00), m=1024, tag='Kyber 512')

    sage: primal_usvp(Kyber512)
    rop: ≈2^140.9, red: ≈2^140.9, δ: 1.004111, β: 382, d: 973, tag: usvp
    
    sage: primal_bdd(Kyber512)
    rop: ≈2^135.8, red: ≈2^134.7, svp: ≈2^135.0, β: 361, η: 406, d: 948, tag: bdd

    sage: params = LWEParameters(n=512, q=3329, Xs=ND.UniformMod(3), Xe=ND.CentredBinomial(eta=8), m=1024)
    sage: primal_usvp(params)
    rop: ≈2^149.0, red: ≈2^149.0, δ: 1.003914, β: 410, d: 944, tag: usvp
         
Status
------

We do not have feature parity with the old estimator yet:

- [X] primal attack on LWE via uSVP using BKZ simulators or assuming the GSA
- [X] primal attack on LWE (effectively) via BDD using BKZ simulators or assuming the GSA
- [ ] hybrid primal attack on LWE
- [ ] dual attack on LWE
- [ ] Aroroa-GB attack on LWE
- [ ] Coded-BKW attack on LWE

We also plan:

- [ ] attack on NTRU pulic keys (using overstretched parameters)
- [ ] SIS attack  
         
Evolution
---------

This code is evolving, new results are added and bugs are fixed. Hence, estimations from earlier
versions might not match current estimations. This is annoying but unavoidable at present. We
recommend to also state the commit that was used when referencing this project.

Contributions
-------------

Our intent is for this estimator to be maintained by the research community. For example, we
encourage algorithm designers to add their own algorithms to this estimator and we are happy to help
with that process.

More generally, all contributions such as bugfixes, documentation and tests are welcome. Please go
ahead and submit your pull requests. Also, don’t forget to add yourself to the list of contributors
below in your pull requests.

At present, this estimator is maintained by Martin Albrecht. Contributors are:

- Benjamin Curtis
- Cedric Lefebvre
- Fernando Virdia
- Florian Göpfert
- James Owen
- Léo Ducas
- Markus Schmidt
- Martin Albrecht
- Rachel Player
- Sam Scott

Citing
------

If you use this estimator in your work, please cite

    | Martin R. Albrecht, Rachel Player and Sam Scott. *On the concrete hardness of Learning with Errors*.
    | Journal of Mathematical Cryptology. Volume 9, Issue 3, Pages 169–203, ISSN (Online) 1862-2984,
    | ISSN (Print) 1862-2976 DOI: 10.1515/jmc-2015-0016, October 2015

A pre-print is available as

    Cryptology ePrint Archive, Report 2015/046, 2015. https://eprint.iacr.org/2015/046

An updated version of the material covered in the above survey is available in
`Rachel Player's PhD thesis <https://pure.royalholloway.ac.uk/portal/files/29983580/2018playerrphd.pdf>`__.

License
-------

The esimator is licensed under the `LGPLv3+ <https://www.gnu.org/licenses/lgpl-3.0.en.html>`__ license.

Acknowledgements
----------------

This project was supported through the European Union PROMETHEUS project (Horizon 2020 Research and
Innovation Program, grant 780701), EPSRC grant EP/P009417/1 and EPSRC grant EP/S020330/1.

Parameters from the Literature
------------------------------

*TODO*
