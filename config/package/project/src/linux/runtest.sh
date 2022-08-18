#!/bin/bash
# Script will run flow123d with given parameters

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RUNTEST=/opt/flow123d/bin/python/runtest.py

# run fterm.sh which will call docker
${SCRIPT_DIR}/../bin/fterm.sh python3 ${RUNTEST} "$@"