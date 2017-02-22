#!/bin/bash
#
# Usage:
# git bisect bisect_test.sh <runtest arguments>
#
# For every bisect level try to build:
# Result in git bisect skip if build fails.
# If build succeeds, call runtest from root directory.
ROOT_DIR="${0%/*}/.."
RUNTEST="${ROOT_DIR}/tests/runtest"

if which pmake; then
    MAKE="pmake"
else
    MAKE="make"
fi

TEST_EXIT_CODE=1
if ${MAKE} all
then    
    if ${RUNTEST} $@
    then 
        TEST_EXIT_CODE=0
    fi    
else
    # force git bisect skip
    exit 125
fi    

exit $TEST_EXIT_CODE