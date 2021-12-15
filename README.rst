Security Estimates for Lattice Problems
=======================================

.. image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Fprompt.ipynb
.. image:: https://readthedocs.org/projects/lattice-estimator/badge/?version=latest
 :target: https://lattice-estimator.readthedocs.io/en/latest/?badge=latest
 :alt: Documentation Status

This `Sage <http://sagemath.org>`__ module provides functions for estimating the concrete security
of
`Learning with Errors <https://en.wikipedia.org/wiki/Learning_with_errors>`__,
`SIS <https://en.wikipedia.org/wiki/Short_integer_solution_problem>`__ and
`NTRU <https://en.wikipedia.org/wiki/NTRU>`__ instances.

The main purpose of this estimator is to give designers an easy way to choose parameters resisting
known attacks and to enable cryptanalysts to compare their results and ideas with other techniques
known in the literature.

Quick Start
-----------

- Usage ::

    >>> from estimator import *
    >>> Kyber512
    LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.00), m=1024, tag='Kyber 512')

    >>> primal_usvp(Kyber512)
    rop: ≈2^141.2, red: ≈2^141.2, δ: 1.004111, β: 382, d: 973, tag: usvp
    
    >>> primal_bdd(Kyber512)
    rop: ≈2^137.8, red: ≈2^136.5, svp: ≈2^137.1, β: 365, η: 400, d: 981, tag: bdd

- `Try it in your browser <https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Fprompt.ipynb>`__.
- `Read the documentation <https://lattice-estimator.readthedocs.io/en/latest/>`__.
  
Status
------

We do not have feature parity with the `old estimator <https://bitbucket.org/malb/lwe-estimator/src/master/>`__ yet:

- ``[x]`` |lwe-primal-binder| :doc:`Primal attack on LWE <../algorithms/lwe-primal>` 
- ``[x]`` |lwe-bkw-binder| :doc:`Coded-BKW attack on LWE <../algorithms/lwe-bkw>` 
- ``[X]`` |gb-binder| :doc:`Aroroa-GB attack on LWE <../algorithms/gb>`
- ``[X]`` :doc:`Dual attack on LWE <../algorithms/lwe-dual>`

.. |lwe-primal-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Flwe-primal.ipynb

.. |lwe-bkw-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Flwe-bkw.ipynb

.. |gb-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/malb/lattice-estimator/jupyter-notebooks?labpath=..%2F..%2Ftree%2Fgb.ipynb
            
We are also planning:

- ``[ ]`` Attacks on NTRU pulic keys (using overstretched parameters)
- ``[ ]`` SIS attacks
         
Evolution
---------

This code is evolving, new results are added and bugs are fixed. Hence, estimations from earlier
versions might not match current estimations. This is annoying but unavoidable. We recommend to also
state the commit that was used when referencing this project.

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
