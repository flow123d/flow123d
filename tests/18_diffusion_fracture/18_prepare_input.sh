#! /bin/bash	
echo
echo "##########      18_diffusion_fracture: MESH and NGH and BCD          ## Start ##"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test18"

#meshing geo file
#gmsh -2 -o input/$NAME.msh -format msh $NAME.geo


$NGH ./ngh.ini
$BCD ./bcd.ini

echo
echo "##########      18_diffusion_fracture: MESH and NGH and BCD          ### End ###"
echo


