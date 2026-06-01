# Running lattice estimator in a Ubuntu 24.04 container

This docker build will run the lattice-estimator using a Ubuntu24.04 LTS base image with the latest stable SageMath installed on Conda. 

## Dependencies

- [Docker](https://docs.docker.com/engine/install/) 

## Build and Run

```bash
docker build . -t lattice-estimator -f docker/ubuntu24.04/Dockerfile
docker run -it lattice-estimator:latest
```

Enter the Python interactive shell to get started running the lattice-estimator:

```console
(lattice_estimator) root@9af6e51a3e65:/opt/lattice_estimator# python
Python 3.12.5 | packaged by conda-forge | (main, Aug  8 2024, 18:36:51) [GCC 12.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from sage.all import *
>>> from estimator import *
>>> schemes.Kyber512
LWEParameters(n=512, q=3329, Xs=D(σ=1.22), Xe=D(σ=1.22), m=512, tag='Kyber 512')
```

## Run pytests

```console
>>> exit()   # exit python interactive shell

(lattice_estimator) root@dd916359f4c6:/opt/lattice_estimator# PYTHONIOENCODING=UTF-8 PYTHONPATH=`pwd` python -m pytest

/opt/conda/envs/lattice_estimator/lib/python3.12/site-packages/nbmake/pytest_plugin.py:6: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
  import pkg_resources
==================================================================== test session starts ====================================================================
platform linux -- Python 3.12.5, pytest-8.4.1, pluggy-1.6.0
rootdir: /opt/lattice_estimator
configfile: pyproject.toml
plugins: anyio-4.10.0, nbmake-1.3.4
collected 85 items                                                                                                                                          

README.rst .                                                                                                                                          [  1%]
docs/schemes/hes.rst .                                                                                                                                [  2%]
docs/schemes/nist-pqc-round-3.rst .                                                                                                                   [  3%]
estimator/cost.py ....                                                                                                                                [  8%]
estimator/gb.py ..                                                                                                                                    [ 10%]
estimator/lwe.py ..                                                                                                                                   [ 12%]
estimator/lwe_bkw.py .                                                                                                                                [ 14%]
estimator/lwe_dual.py .                                                                                                                               [ 15%]
estimator/lwe_guess.py ....                                                                                                                           [ 20%]
estimator/lwe_parameters.py ....                                                                                                                      [ 24%]
estimator/lwe_primal.py ..                                                                                                                            [ 27%]
estimator/nd.py .........................                                                                                                             [ 56%]
estimator/ntru.py ..                                                                                                                                  [ 58%]
estimator/ntru_parameters.py ..                                                                                                                       [ 61%]
estimator/ntru_primal.py ...                                                                                                                          [ 64%]
estimator/prob.py .                                                                                                                                   [ 65%]
estimator/reduction.py ....................                                                                                                           [ 89%]
estimator/simulator.py .                                                                                                                              [ 90%]
estimator/sis.py ..                                                                                                                                   [ 92%]
estimator/sis_lattice.py .                                                                                                                            [ 94%]
estimator/sis_parameters.py .                                                                                                                         [ 95%]
estimator/util.py ..                                                                                                                                  [ 97%]
param_sweep.py ..                                                                                                                                     [100%]

============================================================== 85 passed in 438.12s (0:07:18) ===============================================================
```

