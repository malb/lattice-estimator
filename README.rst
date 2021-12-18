Security Estimates for Lattice Problems
=======================================

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Fprompt.ipynb
.. image:: https://readthedocs.org/projects/lattice-estimator/badge/?version=latest
 :target: https://lattice-estimator.readthedocs.io/en/latest/?badge=latest
 :alt: Documentation Status

This `Sage <http://sagemath.org>`__ module provides functions for estimating the concrete security of `Learning with Errors <https://en.wikipedia.org/wiki/Learning_with_errors>`__ instances.

The main purpose of this estimator is to give designers an easy way to choose parameters resisting known attacks and to enable cryptanalysts to compare their results and ideas with other techniques known in the literature.

Quick Start
-----------

- Usage ::

    >>> from estimator import *
    >>> Kyber512
    LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.00), m=1024, tag='Kyber 512')

    >>> LWE.primal_usvp(Kyber512)
    rop: ≈2^141.2, red: ≈2^141.2, δ: 1.004111, β: 382, d: 973, tag: usvp
    
    >>> r = LWE.estimate.rough(Kyber512)
    usvp                 :: rop: ≈2^111.5, red: ≈2^111.5, δ: 1.004111, β: 382, d: 973, tag: usvp
    dual_hybrid          :: rop: ≈2^134.5, mem: ≈2^130.2, m: 512, red: ≈2^134.3, δ: 1.003611, β: 460, d: 1009, ...

    >>> r = LWE.estimate(Kyber512)
    bkw                  :: rop: ≈2^167.2, m: ≈2^155.1, mem: ≈2^156.1, b: 13, t1: 0, t2: 16, ℓ: 12, #cod: 444, #top: 1, ...
    usvp                 :: rop: ≈2^141.2, red: ≈2^141.2, δ: 1.004111, β: 382, d: 973, tag: usvp
    bdd                  :: rop: ≈2^137.8, red: ≈2^136.8, svp: ≈2^136.8, β: 366, η: 399, d: 978, tag: bdd
    dual                 :: rop: ≈2^165.3, mem: ≈2^127.5, m: 581, red: ≈2^165.2, δ: 1.003573, β: 467, d: 1092, ...
    dual_hybrid          :: rop: ≈2^157.7, mem: ≈2^153.6, m: 512, red: ≈2^157.4, δ: 1.003726, β: 440, d: 1008, ...
    

- `Try it in your browser <https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Fprompt.ipynb>`__.
- `Read the documentation <https://lattice-estimator.readthedocs.io/en/latest/>`__.
  
Status
------

We have feature parity with the `old estimator <https://bitbucket.org/malb/lwe-estimator/src/master/>`__:

- ``[x]`` |lwe-primal-binder| :doc:`primal attacks on LWE <../algorithms/lwe-primal>` 
- ``[X]`` |lwe-dual-binder| :doc:`dual attacks on LWE <../algorithms/lwe-dual>`
- ``[x]`` |lwe-bkw-binder| :doc:`Coded-BKW attack on LWE <../algorithms/lwe-bkw>` 
- ``[X]`` |gb-binder| :doc:`Aroroa-GB attack on LWE <../algorithms/gb>`

.. |lwe-primal-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Flwe-primal.ipynb

.. |lwe-dual-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Flwe-dual.ipynb

.. |lwe-bkw-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Flwe-bkw.ipynb

.. |gb-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Fgb.ipynb
            
but we are also planning:

- ``[ ]`` attacks on `NTRU <https://en.wikipedia.org/wiki/NTRU>`__ pulic keys (using overstretched parameters)
- ``[ ]`` attack on `SIS <https://en.wikipedia.org/wiki/Short_integer_solution_problem>`__ instances
         
Evolution
---------

This code is evolving, new results are added and bugs are fixed. Hence, estimations from earlier
versions might not match current estimations. This is annoying but unavoidable. We recommend to also
state the commit that was used when referencing this project.

.. warning :: We give no API/interface stability guarantees. We try to be mindful but we may reorganize the code without advance warning.

Bugs
----

Please report bugs through the `GitHub issue tracker <https://github.com/malb/lattice-estimator/issues>`__.

Contributions
-------------

At present, this estimator is maintained by Martin Albrecht. Contributors are:

- Benjamin Curtis
- Cedric Lefebvre
- Fernando Virdia
- Florian Göpfert
- James Owen
- Léo Ducas
- Markus Schmidt
- Martin Albrecht
- Michael Walter
- Rachel Player
- Sam Scott

 See :doc:`Contributing <../contributing>` for details on how to contribute.

Citing
------

If you use this estimator in your work, please cite

    | Martin R. Albrecht, Rachel Player and Sam Scott. *On the concrete hardness of Learning with Errors*.
    | Journal of Mathematical Cryptology. Volume 9, Issue 3, Pages 169–203, ISSN (Online) 1862-2984,
    | ISSN (Print) 1862-2976 DOI: 10.1515/jmc-2015-0016, October 2015

A pre-print is available as

    | Cryptology ePrint Archive, Report 2015/046, 2015. https://eprint.iacr.org/2015/046

An updated version of the material covered in the above survey is available in
`Rachel Player's PhD thesis <https://pure.royalholloway.ac.uk/portal/files/29983580/2018playerrphd.pdf>`__.

License
-------

The estimator is licensed under the `LGPLv3+ <https://www.gnu.org/licenses/lgpl-3.0.en.html>`__ license.

Acknowledgements
----------------

This project was supported through the European Union PROMETHEUS project (Horizon 2020 Research and
Innovation Program, grant 780701), EPSRC grant EP/P009417/1 and EPSRC grant EP/S020330/1, and by 
`Zama <https://zama.ai/>`__.
