#! /bin/bash	
echo
echo "##########        01_steady_flow_123d: MESH and NGH and BCD          ## Start ##"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test1"

#meshing geo file
gmsh -2 -algo front2d -3 -algo front3d -o input/$NAME.msh -format msh $NAME.geo

$NGH ./ngh.ini
$BCD ./bcd.ini

echo
echo "##########        01_steady_flow_123d: MESH and NGH and BCD          ### End ###"
echo


