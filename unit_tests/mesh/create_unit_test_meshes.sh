#!/bin/bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
OUTPUT_DIR="${SCRIPTPATH}/meshes"
mkdir ${OUTPUT_DIR}

# target element counts:
# SMALL     3 000
# MEDIUM   30 000
# BIG     300 000

# input geo variants
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

for i in ${!geos[@]}; do
  filename="${geos[$i]%.*}"
  echo "VARIANT: ${filename}"
  # refined small
  mshfile="${OUTPUT_DIR}/${filename}_small_refined.msh"
  gmsh -setnumber fine_step ${steps1_small[$i]} -setnumber mesh ${steps2_small[$i]} -${dim[$i]} -format msh2 -o ${mshfile} ${geos[$i]}
  # refined medium
  mshfile="${OUTPUT_DIR}/${filename}_medium_refined.msh"
  gmsh -setnumber fine_step ${steps1_medium[$i]} -setnumber mesh ${steps2_medium[$i]} -${dim[$i]} -format msh2 -o ${mshfile} ${geos[$i]}
  # refined big
  mshfile="${OUTPUT_DIR}/${filename}_big_refined.msh"
  gmsh -setnumber fine_step ${steps1_big[$i]} -setnumber mesh ${steps2_big[$i]} -${dim[$i]} -format msh2 -o ${mshfile} ${geos[$i]}
  uniform small
  mshfile="${OUTPUT_DIR}/${filename}_small_uniform.msh"
  gmsh -setnumber fine_step ${steps3_small[$i]} -setnumber mesh ${steps3_small[$i]} -${dim[$i]} -format msh2 -o ${mshfile} ${geos[$i]}
  # uniform medium
  mshfile="${OUTPUT_DIR}/${filename}_medium_uniform.msh"
  gmsh -setnumber fine_step ${steps3_medium[$i]} -setnumber mesh ${steps3_medium[$i]} -${dim[$i]} -format msh2 -o ${mshfile} ${geos[$i]}
  # uniform big
  mshfile="${OUTPUT_DIR}/${filename}_big_uniform.msh"
  gmsh -setnumber fine_step ${steps3_big[$i]} -setnumber mesh ${steps3_big[$i]} -${dim[$i]} -format msh2 -o ${mshfile} ${geos[$i]}
  exit 0
done
