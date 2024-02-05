#!/bin/bash
# Usage:
#  
#     auto_build.sh <build_type> <environment> <target_image> <release_tag> [push]
# 

# Stop on first error.
set -e
set -x

build_type=$1; shift            # rel | dbg

environment=$1; shift           # gnu ( | intel)

image_name_base=$1; shift       # flow123d (release build) | ci (continuous integration build, commit hash part of the resutulting image name)
if [ "${image_name_base}" == "flow123d" ];
then
    # Relase image
    target_image="flow123d-${environment}"
else
    # CI image
    target_image=${image_name_base}-${environment}  # ci-${environment}
fi

release_tag=$1; shift                  # version


if [ "$1" == "push" ]
then
    docker_push="docker push"
else
    docker_push="echo docker push"
fi



#######################################
# Paths on host
# TODO: move both version and image_tag to the same location in package
flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}
#destination="`pwd`/publish_${environment}"

################################
# paths inside docker container
flow_install_location=/opt/flow123d
flow_repo_location=`pwd`


#####################
# Docker calls
build_container=build-${build_type}-${environment}




###############################
# Conatiners and images names and tags
imagesversion=`cat ${flow_repo_host}/config/build/image_tag`
release_version=`cat ${flow_repo_host}/version`      

# TODO: pass build type as parameter
#build_image="flow-dev-${environment}-${build_type}:${imagesversion}"








echo "Variables summary"
echo "build_container: '${build_container}'"
echo "release_version: '${release_version}'"
echo "environment: '${environment}'"
echo "imagesversion: '${imagesversion}'"
echo "release_tag: '${release_tag}'"
echo "target_image: '${target_image}'"



######################################################################################################### setup container

# docker rm -f  || echo "container not running"
bin/fterm update
bin/fterm ${build_type}_${environment} --detach ${build_container} 

# we do not know why it later fails at "git config --global --add safe.directory" cmd
# when user is passed
dexec="docker exec -u $(id -u):$(id -g) ${build_container}"      # execute command which will follow
#dexec="docker exec ${build_container}"      # execute command which will follow

#dcp="docker cp ${build_container}"          # Copy files/folders between a container and the local filesystem

# function dexec_setvars_make {
#     # dexec_setvars_make  <make parameters>if [ -s build_tree ]; echo "link";fi
#     ${dexec} /bin/bash -c "cd ${flow_repo_location} && bin/setvars.sh && make $*"
# }


######################################################################################################### build flow123d install container

# copy config
#${dexec} ls ${flow_repo_location}
cp config/config-jenkins-docker-${build_type}.cmake config.cmake

# add full version
echo "set(FLOW_MANUAL_VERSION ${release_tag})" >> config.cmake
${dexec} printenv "HOME"
${dexec} id
${dexec} make -C ${flow_repo_location} set-safe-directory
${dexec} git config --global --add safe.directory '*'

# compile
#DEBUG=--debug=j
${dexec} make -C ${flow_repo_location} clean-all
${dexec} make -C ${flow_repo_location} ${DEBUG} -j4 all
echo "build result: $?"

# Basic test of working binary.
${dexec} bin/flow123d --version
