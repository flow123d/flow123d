#/bin/bash
ABS_FLOW_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )"/.. && pwd )"

# Create benchmark meshes.
# Using docker run with GEOMOP image, call on host.

git_branch=`git rev-parse --abbrev-ref HEAD`
build_dir=build-${git_branch}
${ABS_FLOW_DIR}/config/build/create_gmsh_image.sh
${ABS_FLOW_DIR}/config/build/create_unit_test_meshes.sh ${build_dir}/benchmark_meshes
