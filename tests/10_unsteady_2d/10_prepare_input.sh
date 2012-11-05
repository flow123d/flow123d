#! /bin/bash	
echo
echo "##########         10_unsteady_flow_2d: MESH and NGH and BCD          ## Start ##"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test10"

#generates a little bit different mesh then the original one, 
#check also the parameter in the .geo file
#meshing geo file
#gmsh -2 -o input/$NAME.msh -format msh $NAME.geo

$NGH ./ngh.ini
$BCD ./bcd.ini

echo
echo "##########         10_unsteady_flow_2d: MESH and NGH and BCD          ### End ###"
echo


