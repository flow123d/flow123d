#!/bin/bash
# Script will run flow123d with given parameters

CWD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# run docker shell within current dirrectory
"$CWD/fterm.sh" //opt/flow123d/bin/flow123d "$@"
