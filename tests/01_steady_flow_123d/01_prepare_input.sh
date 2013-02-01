#! /bin/bash	
echo
echo "##########        01_steady_flow_123d: MESH and NGH and BCD          ## Start ##"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test1_new_fbc"

#meshing geo file
#gmsh -2 -algo front2d -3 -algo front3d -o input/$NAME.msh -format msh $NAME.geo

$NGH ./ngh_new_fbc.ini
$BCD ./bcd_new_fbc.ini

echo
echo "##########        01_steady_flow_123d: MESH and NGH and BCD          ### End ###"
echo


