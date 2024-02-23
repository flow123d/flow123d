#!/bin/bash
# Usage:
#  
#     make_packages.sh <environment> <target_image> <release_tag> [push]
# 

# Stop on first error.
set -e
set -x


environment=$1          # gnu
image_name_base=$2      # flow123d | ci
release_tag=$3          # version
if [ "$4" == "push" ]
then
    docker_push="docker push"
else
    docker_push="echo docker push"
fi
# 
#######################################
# Paths on host
# TODO: move both version and image_tag to the same location in package
flow_repo_host="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/../.. && pwd )"
cd ${flow_repo_host}
destination="`pwd`/publish_${environment}"

################################
# paths inside docker container
flow_install_location=/opt/flow123d
flow_repo_location=`pwd`


#####################
# Docker calls
build_type=rel # consistent for fterm and default config files
build_container=contgnurelease




###############################
# Conatiners and images names and tags
imagesversion=`cat ${flow_repo_host}/config/build/image_tag`
#release_version=`cat ${flow_repo_host}/version`      



git_hash=`git rev-parse --short=6 HEAD`
git_branch=`git rev-parse --abbrev-ref HEAD`
build_dir_host=build-${git_branch}


if [ "${image_name_base}" == "flow123d" ];
then
    # Relase image
    target_image="flow123d-${environment}"
else
    # CI image
    target_image=${image_name_base}-${environment}  # ci-${environment}
fi




echo "Variables summary"
echo "build_container: '${build_container}'"
#echo "release_version: '${release_version}'"
echo "environment: '${environment}'"
echo "imagesversion: '${imagesversion}'"
echo "release_tag: '${release_tag}'"
echo "target_image: '${target_image}'"



######################################################################################################### setup container
# Recreate build files and dirs
rm -rf ${build_dir_host}
mkdir ${build_dir_host} && tar xf build_dir.tar -C ${build_dir_host} --strip-components 1
rm -f build_tree
ln -s ${build_dir_host} build_tree
cp ${build_dir_host}/_config.cmake config.cmake


make update-submodules


# docker rm -f  || echo "container not running"
bin/fterm update
bin/fterm rel_$environment --detach ${build_container} #-v `pwd`:${flow_repo_location} 
#docker rm -f ${build_container}
#bin/fterm rel_$environment -- -di --name ${build_container} --volume `pwd`:${flow_repo_location}
#bin/fterm rel_$environment -- --privileged -di --name ${build_container} --volume `pwd`:${flow_repo_location}

dexec="docker exec ${build_container}"      # execute command which will follow
dcp="docker cp ${build_container}"          # Copy files/folders between a container and the local filesystem

function dexec_setvars_make {
    # dexec_setvars_make  <make parameters>if [ -s build_tree ]; echo "link";fi
    ${dexec} /bin/bash -c "cd ${flow_repo_location} && bin/setvars.sh && make $*"
}


######################################################################################################### build flow123d install container

${dexec} echo "HOME=$HOME"
${dexec} make -C ${flow_repo_location} set-safe-directory
${dexec} git config --global --add safe.directory '*'




dexec_setvars_make package

#${dexec} ls ${flow_repo_location}/build_tree/
#${dexec} ls ${flow_repo_location}/build_tree/_CPack_Packages
#${dexec} ls ${flow_repo_location}/build_tree/_CPack_Packages/Linux
${dexec} ls ${flow_repo_location}/build_tree/_CPack_Packages/Linux/TGZ

#${dexec} make -C ${{flow_repo_location}} FORCE_DOC_UPDATE=1 ref-doc
#${dexec} make -C ${{flow_repo_location}} html-doc
#${dexec} make -C ${{flow_repo_location}} doxy-doc

# Local install of source Python packages
${dexec} pip install --user -r ${flow_repo_location}/config/build/requirements.txt

############################################################################################# docker image
install_image="install-${environment}:${imagesversion}"
target_tagged=flow123d/${target_image}:${release_tag}

tmp_install_dir=${flow_repo_location}/build_tree/_CPack_Packages/Linux/TGZ/flow123d-base-${release_tag}-Linux

# have to copy package dir out of the mounted volume ${flow_repo_location}
# TODO: try to use the install target

${dexec} cp -r ${tmp_install_dir} /docker_package
# we only use temporary installation and copy it directly from the build image to the target image
# TODO: need to copy the package out of the vmount, but only for testing

docker rmi -f ${target_tagged}
docker rmi -f flow123d/temporary_build
docker commit ${build_container} flow123d/temporary_build
docker pull "flow123d/${install_image}"

build_date=`date -u +"%Y-%m-%dT%H:%M:%SZ"`

#cp ${destination}/${cmake_package_name} project/src/docker/create/default/${cmake_package_name}
docker build \
     --build-arg base_image=flow123d/${install_image} \
     --build-arg source_image=flow123d/temporary_build \
     --build-arg source_location=/docker_package \
     --build-arg flow_version=${release_tag} \
     --build-arg flow_install_location=${flow_install_location} \
     --build-arg git_hash="${git_hash}" \
     --build-arg build_date="${build_date}" \
     --tag ${target_tagged} \
     ${flow_repo_host}/config/build/dockerfile

# Push only if $3 == push
${docker_push} ${target_tagged}

##########################################################################################################šš


# name of the archives from CMake CPack tool
cmake_package_name=Flow123d-${release_tag}-Linux.tar.gz

# final archive names
base_name=flow123d_${release_tag}
#docker_arch_name=${base_name}_docker_image.tar.gz
lin_arch_name=${base_name}_linux_install.tar.gz
win_arch_name=${base_name}_windows_install.exe

# current date in two forms
current_date=`date +"%d-%m-%Y %T"`

# path to the generated pdf
pdf_location=doc/reference_manual/flow123d_doc.pdf
ist_location=doc/reference_manual/input_reference.json




############################################################################################ copy testsmake doc

# create destination folders
mkdir -p ${destination}/tests
# clean and copy tests folder

dexec_setvars_make clean-tests
${dcp}:${flow_repo_location}/tests/. ${destination}/tests

# delete runtest because we will have to create other runtest for docker
rm -rf ${destination}/tests/runtest
############################################################################################ make doc

dexec_setvars_make -j4 all                      # compile flow123d
dexec_setvars_make FORCE_DOC_UPDATE=1 ref-doc   # generate latex doc
dexec_setvars_make html-doc                     # generate html doc
dexec_setvars_make doxy-doc                     # generate source doc

mkdir -p ${destination}/htmldoc
mkdir -p ${destination}/doxygen
#mkdir -p ${destination}/config/docker/
mkdir -p ${destination}/bin/

${dcp}:${flow_repo_location}/build_tree/${pdf_location}             ${destination}/flow123d_${release_tag}_doc.pdf
${dcp}:${flow_repo_location}/build_tree/htmldoc/html/src/.          ${destination}/htmldoc
${dcp}:${flow_repo_location}/build_tree/doc/online-doc/flow123d/.   ${destination}/doxygen
${dcp}:${flow_repo_location}/${ist_location}                        ${destination}/input_reference.json
${dcp}:${flow_repo_location}/bin/fterm                              ${destination}/bin/fterm

echo "${release_tag}" > ${destination}/version
echo "${target_tagged}" > ${destination}/imagename

############################################################################################## Linux package

# TODO: simplify creation of the linux package to single command and run it inside the docker
# ... does not need make, cmake and other tools out of the docker


mkdir -p install-linux
#  mkdir -p make parent directories as needed
cmake \
    -DFLOW_VERSION="${release_tag}" \
    -DFLOW123D_ROOT="${flow_repo_location}" \
    -DIMAGE_TAG="${target_tagged}" \
    -DDEST="${destination}" \
    -S ${flow_repo_host}/config/build/project -B install-linux

make -C install-linux package
mv install-linux/${base_name}.tar.gz ${destination}/${lin_arch_name}
echo "{\"build\": \"${current_date}\", \"hash\": \"${git_hash}\"}" > ${destination}/flow123d_${release_tag}_linux_install.json

############################################################################################## Windows package

mkdir -p install-win
cp -r ${flow_repo_host}/config/build/project/src/windows/* install-win
cp -r ${flow_repo_host}/tests install-win
${dcp}:${flow_repo_location}/build_tree/htmldoc/html/src/.  install-win/htmldoc
${dcp}:${flow_repo_location}/build_tree/${pdf_location}     install-win/flow123d_${release_tag}_doc.pdf


echo "${release_tag}" > install-win/version
echo "${target_tagged}" > install-win/imagename


# The 'docker run' command first creates a writeable container layer over the specified image, and then starts it using the specified command
docker run -i --rm --user ${uid}:${gid} -v `pwd`/install-win:/nsis-project hp41/nsis /nsis-project/install.nsi

mv install-win/${win_arch_name} ${destination}/${win_arch_name}
echo "{\"build\": \"${current_date}\", \"hash\": \"${git_hash}\"}" > ${destination}/flow123d_${release_tag}_windows_install.json
	

        

# list build packages
ls -l ${destination}

############################################################################################## Windows package

