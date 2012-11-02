#! /bin/bash	
echo
echo "##########        13_edge solute mixing: MESH and NGH and BCD          ## Start ##"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test13"

#generates a little bit different mesh then the original one
#meshing geo file
#gmsh -2 -o input/$NAME.msh -format msh $NAME.geo

$NGH ./ngh.ini
$BCD ./bcd.ini

echo
echo "##########        13_edge solute mixing: MESH and NGH and BCD          ### End ###"
echo


