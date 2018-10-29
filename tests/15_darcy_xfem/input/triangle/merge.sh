#!/bin/bash
#
# Compute all tests:
#   run_scale_test.sh 1 
#
# Just process data:
#   run_scale_test.sh
set -x
PWD=`pwd`

gmsh_merge_file="merge_gmsh.geo"

for index in  0 1 2 3 4 5
do
      mesh_1="triangle_${index}.msh"
      mesh_2="well.msh"
      
      output_name="triangle_well_${index}.msh"

      echo "merging: $mesh_1 + $mesh_2"

      printf "Merge \"$mesh_1\";
Merge \"$mesh_2\";
Save \"$output_name\";
Exit;" > $gmsh_merge_file
     gmsh $gmsh_merge_file -
done
