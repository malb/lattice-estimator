FROM sagemath/sagemath:10.3
COPY --chown=sage:sage . ${HOME}
USER sage
RUN /home/sage/sage/local/var/lib/sage/venv-python3.11.1/bin/python3 -m pip install --no-cache-dir notebook jupyterlab
