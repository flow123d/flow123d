#!/bin/bash
# Usage:
# ./run_tests_with_it_image.sh <release_tag> <command> [args]

# Stop on first error
set -e
set -x

environment=$1; shift   # gnu | intel
release_tag=$1; shift   # version	(4.0.3dev_cdd097)
command_with_args=$@    # run the command with arguments

target_image="stepanmoc/it-${environment}:${release_tag}"

if [[ "$(docker images -q ${target_image} 2> /dev/null)" == "" ]]; then
    echo "Docker image '${target_image}' not found. Please build the image first."
    exit 1
fi

container_id=$(docker run -d ${target_image})
docker exec ${container_id} bash -c "cd /opt/flow123d/bin"

for test_dir in "${test_dirs[@]}"; do
  docker exec ${container_id} ${command_with_args}
done

# Wait for tests to finish
wait

# Stop and remove the container
docker stop ${container_id}
docker rm ${container_id}