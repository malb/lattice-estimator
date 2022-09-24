Contributing
============

Our intent is for this estimator to be maintained by the research community. For example, we
encourage algorithm designers to add their own algorithms to this estimator and we are happy to help
with that process.

More generally, all contributions such as bugfixes, documentation and tests are welcome. Please go
ahead and submit your pull requests.

Design & Tips
-------------

We use Python classes as namespaces, something like this::

    class SuperCoolAttack:

        @staticmethod
        def some_complex_function(a, b, c, log_level=5):
            cost = Cost(rop=a*b*c, a=a)
            cost.register_impermanent(rop=True, a=False) # (6)
            return cost
          
        def __call__(self, params, a=None, log_level=5):
            """
            Estimate cost of super cool attack [C:Author21]_. # (8)
            
            :param params: LWE parameters.
            :param a: The famous `a` parameter.              
            
            """
            params = params.normalize() # (1)
            with local_minimum(params.n//2, stop=params.n, log_level=log_level+1) as it:
                for n in it: # (2)
                    cost = self.some_complex_function(a, b, c, log_level=log_level)
                    it.update(cost) # (2): submit current value to pick next one
                cost = it.y # extract final solution
                Logging.log("sca", log_level, f"bla: {repr(cost)}") # (3)
                    
            cost["tag"] = "super-cool-attack" # (4)
            cost["problem"] = params # (5)
            return cost

        __name__ = "super_cool_attack"
            
    super_cool_attack = SuperCoolAttack() # (7)
        
We explain what is going on above:
    
1. LWE objects know how to normalize themselves by calling ``params.normalize()``. We assume that high-level functions (such as ``__call___`` above) call ``params.normalize()``.

2. Often optimizing parameters means finding the optimimum in some range. We provide some syntactical sugar to make this easy/readable.

3. Logging takes levels and modules, for the latter to work you need to add an entry to ``Logging.loggers`` 

4. Costs should indicate where they come from.

5. We assume "costs" know their "problems".

6. When constructing cost objects we should indicate which entries scale with repetitions and which do not. For example, the block size of lattice reduction would not, the overall cost would.

7. The user visible object is then ``super_cool_attack``.

8. Add references to :doc:`References <../references>`.

.. note :: Don’t forget to add yourself to the list of contributors in the README in your pull request.
   
Tests
-----

- The easiest way to run tests is to run ``pytest``, ``python -m pytest`` (for those running inside ``conda`` or ``sage -sh``), or ``sage --python -m pytest`` in the root directory of this repository. These tests are also run on each commit. Code that isn't tested is broken, so ideally every function should have an example doctest.
  - To install the depdencies for testing, run ``pip install -r requirements.txt``.
- We also enforce a coding style using `flake8 <https://flake8.pycqa.org/en/latest/>`__.
- You should also test building the documentation locally before creating a pull request (see below).

Documentation
-------------

The documentation for the ``estimator`` is available `online <https://lattice-estimator.readthedocs.io/>`__ and it can be generated locally by running the following command in the root directory of this repository::

    make html

You will need `sphinx-book-theme <https://sphinx-book-theme.readthedocs.io/en/latest/>`__ installed for this to work. If the documentation was previously generated locally use the following command to regenerate everything from scratch::
    
    make clean && make html

You can produce Jupyter notebooks, run::

    make jupyter

You will need `sphinxcontrib-jupyter <https://github.com/QuantEcon/sphinxcontrib-jupyter>`__ installed for this to work.

Note that `sphinxcontrib-jupyter <https://github.com/QuantEcon/sphinxcontrib-jupyter>`__ currently has a `bug <https://github.com/QuantEcon/sphinxcontrib-jupyter/issues/339>`__ which means ``make html && make jupyter`` will fail. However, ``make html && make clean && make jupyter`` should succeed.
