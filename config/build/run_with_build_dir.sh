# Usage: run_with_build_dir.sh [-c] <environment> <command> [args]
#
# environment= dbg_gnu | rel_gnu | dbg_intel ....
#
# Assumes cloned repository set to correct branch. 
# Extract the build_dir tar ball.
# Run the <command> within the build container with the extracted build dir.
#
# -c .. Continue with curent build_dir, do not extract the tarball.

set -e
set -x

if [ "$1" == "-c" ];then
    PRESERVE_BUILD_DIR=1
    shift
fi

env=$1;shift
case $env in
    dbg_gnu|rel_gnu|dbg_intel|rel_intel):
    ;;
    *)
    echo Wrong environment.
    exit 1
esac


command_with_args=$@

git_branch=`git rev-parse --abbrev-ref HEAD`
build_dir_host=build-${git_branch}

# Recreate build files and dirs
if [ -z "$PRESERVE_BUILD_DIR" ]; then
    rm -rf ${build_dir_host}
    mkdir ${build_dir_host} && tar xf build_dir.tar -C ${build_dir_host} --strip-components 1
    rm -f build_tree
    ln -s ${build_dir_host} build_tree
    cp --preserve=timestamps ${build_dir_host}/_config.cmake config.cmake
fi
make update-submodules
#bin/fterm ${env} --no-term  exec make --dry-run fast-flow123d
    
# run the command
bin/fterm ${env} --no-term  exec ${command_with_args}



# docker ps
# docker restart ${build_container}
# 
# dexec="docker exec ${build_container}"      # execute command which will follow





################################################################################################## build flow123d install container

# 
# copy config
# ${dexec} ls ${flow_repo_location}
# ${dexec} cp ${flow_repo_location}/config/config-jenkins-docker-${config_name}.cmake ${flow_repo_location}/config.cmake
# 
# overload git protection
# set_safe_dir="${dexec} git config --global --add safe.directory"
# for d in \
#     ${flow_repo_location} \
#     ${flow_repo_location}/bin/yaml_converter \
#     ${flow_repo_location}/src/dealii \
#     ${flow_repo_location}/third_party/bparser \
#     ${flow_repo_location}/third_party/json-3.10.5 \
#     ${flow_repo_location}/third_party/gtest-1.10.0
# do
#     ${set_safe_dir} $d
# done
# 
# compile
# ${dexec} make -C ${flow_repo_location} -j4 all
# 
# tar build tree
# ${dexec} cp ${flow_repo_location}/config/config-jenkins-docker-${build_type}.cmake ${flow_repo_location}/build_tree/_config.cmake
# tar -cvf buid_tree.tar build_tree
