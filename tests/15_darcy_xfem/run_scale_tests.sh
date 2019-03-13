#!/bin/bash
#
# Compute all tests:
#   run_scale_test.sh 1 
#
# Just process data:
#   run_scale_test.sh
set -x
PWD=`pwd`


flow123d="/home/paulie/docker_flow123d/flow123d/bin/flow123d"
# flow123d="$HOME/docker_flow123d/flow123d/bin/flow123d"

# CON file 
# template_file="triangle_const_template.yaml"
# template_file="triangle_sigma_template.yaml"
# template_file="sigma_a_template.yaml"
# template_file="sigma_source_a_template.yaml"
# template_file="sigma_2w_template.yaml"
# template_file="sigma_5w_template.yaml"
# template_file="triangle_sigma_source_template.yaml"
# template_file="triangle_sigma_2w_koeppl_template.yaml"
# template_file="tetrahedron3t_sigma_template.yaml"

# template_file="23_2d_circle_well_template.yaml"
# template_file="24_2d_circle_well_source_template.yaml"
# template_file="24_2d_circle_well_source_only_template.yaml"
# template_file="24a_2d_circle_well_source_only_template.yaml"
# template_file="24b_2d_circle_well_source_only_template.yaml"
# template_file="24c_2d_circle_well_source_only_template.yaml"

# template_file="24d_2d_circle_well_source_only_template.yaml"
# template_file="24d_2d_circle_well_source_template.yaml"
# template_file="24d_2d_circle_well_no_source_template.yaml"

# template_file="26_square_2w_template.yaml"
# template_file="26_square_2w_source_template.yaml"
template_file="26_square_2w_source_only_template.yaml"

# template_file="51_single_aquifer_analytical_source_template.yaml"
# template_file="51_single_aquifer_analytical_template.yaml"

# template_file="25_3d_cylinder_well_template.yaml"
# template_file="25_3d_cylinder_well_source_template.yaml"
# template_file="25_3d_cylinder_well_source_only_template.yaml"

# template_file="31_block_2w_template.yaml"
# template_file="31_block_2w_source_template.yaml"
# template_file="31_block_2w_source_only_template.yaml"

# template_file="31_block_5w_template.yaml"
# template_file="31_block_5w_source_template.yaml"

# template_file="33_block_2w_template.yaml"
# template_file="33_block_2w_source_template.yaml"

yaml_file="scale_test_run.yaml"

summary_file="summary_${template_file}.log"

touch $summary_file
echo " index p2D_diff v2D_diff p3D_diff v3D_diff enr_v_dof sing_lp meshfile"  >> $summary_file
for index in  1 2 3 4 5 6
do
  
#       mesh_name="triangle_const_${index}.msh"
#       mesh_name="big_circle_${index}.msh"
#     mesh_name="triangle\/triangle_well_${index}.msh"

#     mesh_name="big_circle_conv\/big_circle_well_${index}.msh"
#     mesh_name="big_circle_conv_well\/big_circle_well_${index}.msh"
#     mesh_name="big_circle_conv\/big_circle_${index}.msh"

    mesh_name="square\/square_${index}.msh"
#     mesh_name="square\/square_2w_${index}.msh"
#     mesh_name="square\/square_5w_${index}.msh"

#     mesh_name="cylinder_conv\/cylinder_well_${index}.msh"
#     mesh_name="cylinder_conv\/cylinder_${index}.msh"

# mesh_name="block_conv\/block_${index}.msh"
# mesh_name="block_conv\/block_2w_${index}.msh"
# mesh_name="block_conv\/block_5w_${index}.msh"

# mesh_name="block_conv_closer\/block_2w_${index}.msh"
    
#     mesh_name="aquifers\/single_aquifer_12d_${index}.msh"
#       mesh_name="aquifer_${index}.msh"
#       mesh_name="cylinder_${index}.msh"
#       mesh_name="tetrahedron3t_${index}.msh"
      echo $mesh_name
      output_name="output_${index}"
      
      max_level=`expr 6 - $index`
      if ((max_level<1)); then max_level=1; fi
      
#       max_level=1

      cat $template_file | sed "s/MESH_FILE_PARAM/${mesh_name}/; s/MAX_LEVEL_PARAM/${max_level}/" > $yaml_file
      #cat $yaml_file | sed "s/MAX_LEVEL_PARAM/${max_level}/" > $yaml_file
      
      time $flow123d -s $yaml_file --no_profiler -i . -o "${output_name}" > flow123d.out
      
      p_2D_err=`cat "${output_name}/solution_error" | grep 'pressure error 2d:' | sed "s/.*://"`
      v_2D_err=`cat "${output_name}/solution_error" | grep 'velocity error 2d:' | sed "s/.*://"`
      p_3D_err=`cat "${output_name}/solution_error" | grep 'pressure error 3d:' | sed "s/.*://"`
      v_3D_err=`cat "${output_name}/solution_error" | grep 'velocity error 3d:' | sed "s/.*://"`
      enr_v_dof=`cat "${output_name}/solution_error" | grep 'enr vel dof:' | sed "s/.*://"`
      sing_lp=`cat "${output_name}/solution_error" | grep 'sing LP:' | sed "s/.*://"`
      echo "${index}   ${p_2D_err}   ${v_2D_err}   ${p_3D_err}   ${v_3D_err}  ${enr_v_dof}   ${sing_lp}   ${mesh_name}" >> $summary_file
done
