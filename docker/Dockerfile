# A Dockerfile that can be used in isolation
# syntax=docker/dockerfile:1
FROM sagemath/sagemath:latest

user root
RUN apt update && apt upgrade -y
RUN apt install -y git

RUN git clone https://github.com/malb/lattice-estimator.git /lattice-estimator
WORKDIR /lattice-estimator

RUN sage -pip install -r /lattice-estimator/requirements.txt
CMD ["sage"]
