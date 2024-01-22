#!/bin/bash
set -x
# Prepare benchmark meshes using GMSH in a docker image `geomop-gnu`. 
# Syntax:
# 
#    create_unit_test_meshes.sh <OUTPUT_DIR> [<MAX_SIZE>]
#
#  Create meshes in the directory <OUTPUT_DIR>. 
#  MAX_SIZE is one of [small|medium|big], with 'big' default. 


# suppose that
# - pwd is the root dir of flow123d repository
# - OUTPUT_DIR is also under the root dir of flow123d repository

# SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
# cd ${SCRIPTPATH}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
DEFAULT_OUTPUT_DIR="${SCRIPT_DIR}/../../build_tree/benchmark_meshes"
OUTPUT_DIR="${1:-${DEFAULT_OUTPUT_DIR}}"
if [ -z "$OUTPUT_DIR" ]; then
    echo "Missing output dir."
    exit 1
fi

MAX_SIZE=${2:-big}
SIZE_NAMES=("small" "medium" "big")
for i_size in ${!SIZE_NAMES[@]}; do
    if [ "${MAX_SIZE}" == "${SIZE_NAMES[$i_size]}" ]; then
        MAX_SIZE_INT=$i_size
    fi
done
    
mkdir -p ${OUTPUT_DIR}
echo "pwd: $(pwd)"

# gmsh in Docker container
docker_container="flow123d/geomop-gnu:2.0.0"
docker_contname="geomop-gnu"
echo "start docker container: '${docker_container}'"
docker run -t -d --name ${docker_contname} -w /$(pwd) -v /$(pwd):/$(pwd) ${docker_container}
dexec="docker exec -u $(id -u):$(id -g) ${docker_contname}"
gmsh="${dexec} gmsh -format msh2"

# target element counts:
# SMALL     3 000
# MEDIUM   30 000
# BIG     300 000

# input geo variants
ut_dir="unit_tests/coupling"
geos=("lshape_2D.geo" \
      "lshape_3D.geo" \
      "square_2D.geo" \
      "cube_3D.geo")

# corresponding dimension
dim=(2 3 2 3)


# SMALL fine_step
steps1_small=(1e-6 3.5e-2 1e-6 3.2e-2)
# SMALL mesh
steps2_small=(0.115 0.55 0.155 0.45)
# SMALL (uniform) mesh
steps3_small=(0.05 0.222 0.028 0.1245)

# MEDIUM fine_step
steps1_medium=(2.2e-8 3.3e-3 4e-6 2.8e-3)
# MEDIUM mesh
steps2_medium=(0.0335 0.6 0.04 0.46)
# MEDIUM (uniform) mesh
steps3_medium=(0.0152 0.097 0.0088 0.054)

# BIG fine_step
steps1_big=(4e-9 2.9e-4 6.03e-8 2.68e-4)
# BIG mesh
steps2_big=(0.0097 0.5 0.014 0.37)
# BIG (uniform) mesh
steps3_big=(0.0048 0.045 0.00275 0.0245)

function make_mesh {
  fine_step=$1
  mesh_step=$2
  dim=$3
  size_int=$4
  if [ "$size_int" -le "$MAX_SIZE_INT" ]; then
      ${gmsh} -setnumber fine_step $fine_step -setnumber mesh ${mesh_step} -${dim[$i]} -o ${mshfile} ${geofile}
  fi
}

function make_mesh_variants {
  i_mesh=$1

  geofile="${ut_dir}/${geos[$i_mesh]}"
  filename="${geos[$i_mesh]%.*}"
  echo "generate meshes for: ${filename}"
  for i_size in 0 1 2; do
    size_str=${SIZE_NAMES[$i_size]}
        
    
    eval refined_step_min=\${steps1_${size_str}[$i_mesh]}
    eval refined_step_max=\${steps2_${size_str}[$i_mesh]}
    eval uniform_step=\${steps3_${size_str}[$i_mesh]}
    
    # refined
    mshfile="${OUTPUT_DIR}/${filename}_refined_${size_str}.msh"
    make_mesh ${refined_step_min} ${refined_step_max} ${dim[$i_mesh]} ${i_size}
    
    # uniform
    mshfile="${OUTPUT_DIR}/${filename}_uniform_${size_str}.msh"
    make_mesh ${uniform_step} ${uniform_step} ${dim[$i_mesh]} ${i_size}
  done
}

for i in ${!geos[@]}; do
  make_mesh_variants $i
done

docker stop ${docker_contname}
docker rm ${docker_contname}
