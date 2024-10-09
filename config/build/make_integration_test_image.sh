#!/bin/bash
# Usage:
# ./make_integration_test_image.sh <environment> <image_source_name_base> <release_tag> [push]

# Stop on first error
set -e
set -x

environment=$1              # gnu
image_source_name_base=$2   # flow123d | ci
release_tag=$3              # version
docker_push=$4              # push

docker_push="docker push"

# Paths
flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}

# Docker calls
build_container=contgnurelease
source_image="stepanmoc/${image_source_name_base}-${environment}:${release_tag}"
target_image="stepanmoc/it-${environment}:${release_tag}"

docker rmi -f ${target_tagged}
docker pull ${source_image}

# Build a new Docker image by adding VTK and copying the test directory
docker build --build-arg source_image=${source_image} \
             --build-arg flow_version=${release_tag} \
             --build-arg git_hash=$(git rev-parse --short HEAD) \
             --build-arg build_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ") \
             --tag ${target_image} \
             ${flow_repo_host}/config/build/it_docker_file

# Push the image to Docker Hub
docker push ${target_image}

# Done
echo "Docker image '${target_image}' has been built and pushed."
