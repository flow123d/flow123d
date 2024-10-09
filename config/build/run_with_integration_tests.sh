#!/bin/bash
# Usage: run_in_it_gnu.sh [-c] <environment> <command> [args...]
#
# environment = dbg_gnu | rel_gnu | dbg_intel | rel_intel ...
# <command> and [args] are the command and its arguments to run inside the it-gnu image
#
# -c .. Continue with the current build_dir, do not extract the tarball.

set -e
set -x

env=$1; shift
command_with_args=$@

# Set paths and variables
flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
image_tag="stepanmoc/it-gnu:${env}"

docker pull ${image_tag}

docker run --rm \
    -v ${flow_repo_host}:/flow123d \
    -w /flow123d \
    ${image_tag} \
    ./config/build/run_with_build_dir.sh ${env} ${command_with_args}
