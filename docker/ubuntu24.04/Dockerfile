FROM ubuntu:24.04

# Disabling installing extra packages with apt-get install
RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker

RUN DEBIAN_FRONTEND=noninteractive apt-get update && \
  apt-get install -y \
  make \
  wget \
  python3 \
  python3-pip \
  && rm -rf /var/lib/apt/lists/*
  
# download conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

# Add conda to PATH
ENV PATH=/opt/conda/bin:$PATH

# Run conda init
RUN conda init

# Restart shell
RUN bash -l

# Install sage 
# add conda-forge channel
RUN conda config --add channels conda-forge 
RUN conda config --set channel_priority strict

# accept terms of service
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

# create conda environment with sage and python    
RUN conda create -y -n lattice_estimator sage python=3.12.5
RUN conda run -n lattice_estimator pip install sphinxcontrib-jupyter

# copy lattice-estimator repo to container
WORKDIR /opt/lattice_estimator
COPY . .

# ensure interactive shells auto-activate the env, activate the conda env, and install the sagemath kernel
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc && \
    echo "conda activate lattice_estimator" >> /root/.bashrc && \
    /opt/conda/envs/lattice_estimator/bin/sage -python -m ipykernel install --name sagemath --display-name "SageMath" --sys-prefix

# Install pytest dependencies
RUN conda install -n lattice_estimator -c conda-forge numpy nbmake pytest traitlets
# dependency not available in conda-forge hence use pip
RUN pip install --upgrade sphinxcontrib-jupyter 
RUN make jupyter
