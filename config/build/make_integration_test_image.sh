#!/bin/bash
# Usage:
# make_integration_test_image.sh <environment> <target_image> <release_tag> [push]


#######################################
# GET INPUT ARGMUNETS
#######################################

# Stop on first error.
set -e
set -x

build_type=$1;         # rel | dbg
environment=$2;        # gnu ( | intel)
image_name_base=$3;    # ci (continuous integration build, commit hash part of the resutulting image name)
release_tag=$4;        # version
can_push=$5;           # push

# TODO: use push or echo check by if [ "$can_push" == "push" ]
# Use push for testing purposes
docker_push="docker push"

# what todo next we can add tests file to docker for testing, but it will be only
# for the local testing. Need to use ci_gnu sa a base image for the testing and 
# add extra libraries for the testing. Like vtk.
# Base image (ci_gnu) should contain correct runtest from _CPack_Packages


#######################################
# SET VARIABLES
#######################################

# Paths on host
flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}
destination="`pwd`/publish_${environment}"

# Paths inside docker container
flow_install_location=/opt/flow123d
flow_repo_location=`pwd`

# Docker calls
build_container=itcontgnurelease

# Image version
imagesversion=`cat ${flow_repo_host}/config/build/image_tag`

# Git hash and branch
git_hash=`git rev-parse --short=6 HEAD`
git_branch=`git rev-parse --abbrev-ref HEAD`

# Build directory
build_dir_host=build-${git_branch}

# Base image
base_image=${image_name_base}-${environment}
base_tagged=${base_image}:${imagesversion}

# Target image	
target_image=it-${environment}
target_tagged=${target_image}:${imagesversion}



#######################################
# VARIABLES SUMMARY CHECK
#######################################
echo "VARIABLES SUMMARY CHECK"
echo "build_type: ${build_type}"
echo "environment: ${environment}"
echo "image_name_base: ${image_name_base}"
echo "release_tag: ${release_tag}"
echo "can_push: ${can_push}"
echo "docker_push: ${docker_push}"
echo "imagesversion: ${imagesversion}"
echo "git_hash: ${git_hash}"
 
