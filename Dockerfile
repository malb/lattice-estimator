FROM sagemath/sagemath:latest
COPY --chown=sage:sage . ${HOME}
