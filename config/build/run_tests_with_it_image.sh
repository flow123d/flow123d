#!/bin/bash
# Usage:
# ./run_tests_with_it_image.sh <release_tag> <command> [args]

# Stop on first error
set -e
set -x

environment=$1; shift   # gnu | intel
release_tag=$1; shift   # version	(4.0.3dev_cdd097)
command_with_args=$@    # run the command with arguments


flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}

target_image="flow123d/it-${environment}:${release_tag}"

container_id=$(docker run -d --tty=true --interactive=false -v ${flow_repo_host}/tests:/opt/flow123d/bin/tests ${target_image})

trap 'docker stop ${container_id}; docker rm ${container_id}' EXIT

docker exec ${container_id} bash -c "cd /opt/flow123d/bin && ${command_with_args}"