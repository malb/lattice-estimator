FROM sagemathinc/cocalc:latest

SHELL ["/bin/bash", "-c"]
RUN SAGE_ROOT=`pwd`/sage && \
    export SAGE_ROOT="$SAGE_ROOT" && \
    source "$SAGE_ROOT/local/bin/sage-env" && \
    pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache notebook jupyterlab

COPY . ${HOME}
