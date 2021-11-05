FROM sagemath/sagemath:latest

SHELL ["/bin/bash", "-c"]
RUN SAGE_ROOT=`pwd`/sage && \
    export SAGE_ROOT="$SAGE_ROOT" && \
    source "$SAGE_ROOT/local/bin/sage-env" && \
    pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache notebook jupyterlab

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

USER root
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}

COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
