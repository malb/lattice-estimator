# Running lattice estimator in a container

Requires `docker` or `podman` to run the containers.

From the base of the repository, run the following to build an image and run the
tests within the container.

```bash
docker build -t lattice-estimator -f docker/Dockerfile.dev .
docker run -dit --name lattice-estimator-tests lattice-estimator
docker exec lattice-estimator-tests sage -sh -c pytest
```

Note that due to [this open
ticket](https://trac.sagemath.org/ticket/34242#comment:20) on Sage, the
published sage container is using an OEL'ed version of Ubunutu, and so we have
to wait for the sagemath image to be updated in order to use `Dockerfile`
standlone (`git clone`ing from inside the container).

