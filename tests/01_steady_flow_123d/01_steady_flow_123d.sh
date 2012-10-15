#! /bin/bash	
# Use the bash shell to run the NGH sample
echo
echo "################################################################################"
echo "##########                     01_steady_flow_132d                   ## Start ##"
echo "################################################################################"
echo

#choose one of following:	
#	flow
#	flow_mumps 
#	flow_vtk 
#	flow_matis
CONFIG="flow"

NAME="test1"

CON=../../bin/ini2json.sh
NGH=../../bin/ngh/bin/ngh 
BCD=../../bin/bcd/bin/bcd
FLOW=../../bin/flow123d


#meshing geo file
gmsh -2 -algo front2d -3 -algo front3d -o input/$NAME.msh -format msh $NAME.geo

$CON  -ini ./$CONFIG.ini -final ./$CONFIG.con
$NGH ./ngh.ini
$BCD ./bcd.ini
$FLOW  -i input -o output -s ./$CONFIG.con #-ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor


echo
echo "################################################################################"
echo "##########                     01_steady_flow_132d                   ##  End  ##"
echo "################################################################################"
echo


