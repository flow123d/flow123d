#!/bin/bash

#######################################
# Paths on host
flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}


#######################################
# Build image
docker build ${flow_repo_host}/config/build/gmsh_dockerfile -t flow_gmsh
