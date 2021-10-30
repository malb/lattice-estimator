#!/usr/bin/env bash
###############################################################################
#                              Run Sage doctests
###############################################################################
if [ "x$SAGE_ROOT" == "x" ]
then
    SAGE_ROOT=$(sage -c "import os; print(os.environ['SAGE_ROOT'])")
    export SAGE_ROOT="$SAGE_ROOT"
    export PATH="$SAGE_ROOT/bin:$SAGE_ROOT/src/bin:$PATH"
    export SAGE_LOCAL="$SAGE_ROOT/local"
    export DOT_SAGE=$HOME/.sage/
fi

RESULT=0

for file in "$@"; do
    PYTHONIOENCODING=UTF-8 PYTHONPATH=$(pwd) sage-runtests "$file"
    RESULT=$(( RESULT + $? ))
done

exit $RESULT
