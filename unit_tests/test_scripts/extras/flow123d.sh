#!/bin/bash

# get absolute dir in which is this script stored
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python ${SCRIPT_DIR}/flow123d.py "$@"