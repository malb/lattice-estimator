#!/usr/bin/env bash
###############################################################################
#                               Fix-up Doctests
#
# Please don't just blindly call this to make failures go away,
# but review all changes.
###############################################################################
if [ "x$SAGE_ROOT" == "x" ]
then
    SAGE_ROOT=$(sage -c "import os; var = 'SAGE_ROOT' if 'SAGE_ROOT' in os.environ else 'SAGE_VENV'; print(os.environ[var])")
    export SAGE_ROOT="$SAGE_ROOT"
    export PATH="$SAGE_ROOT/bin:$SAGE_ROOT/src/bin:$PATH"
    export SAGE_LOCAL="$SAGE_ROOT/local"
    export DOT_SAGE=$HOME/.sage/
fi

PYTHONIOENCODING=UTF-8 PYTHONPATH=$(pwd) sage-fixdoctests $@
