#!/bin/bash

set -x

    VALGRIND_ARGS=""
    while [ "$1" != "--" ]
    do
        VALGRIND_ARGS="$VALGRIND_ARGS $1"
        shift
    done        
    shift

# call flow123d under valgrind
valgrind $VALGRIND_ARGS --suppressions=${0%/*}/python.supp ${0%/*}/flow123d "$@"