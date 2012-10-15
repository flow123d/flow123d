#! /bin/bash	
# Use the bash shell to run the NGH sample
echo
echo "################################################################################"
echo "##########                     01_steady_flow_132d                   ## Start ##"
echo "################################################################################"
echo

#choose one of following:
#	all		
#	flow
#	flow_vtk 
#	flow_mumps 
#	flow_matis
CONFIG="all"

NAME="test1"

CON=../../bin/ini2json.sh
NGH=../../bin/ngh/bin/ngh 
BCD=../../bin/bcd/bin/bcd
FLOW=../../bin/flow123d


#meshing geo file
gmsh -2 -algo front2d -3 -algo front3d -o input/$NAME.msh -format msh $NAME.geo

$NGH ./ngh.ini
$BCD ./bcd.ini

if [ "$CONFIG" == "all" ]; then
	$CON  -ini ./flow.ini -final ./flow.con
	$CON  -ini ./flow_vtk.ini -final ./flow_vtk.con
	$CON  -ini ./flow_mumps.ini -final ./flow_mumps.con
	$CON  -ini ./flow_matis.ini -final ./flow_matis.con
else
	$CON  -ini ./$CONFIG.ini -final ./$CONFIG.con
fi

if [ $CONFIG == "all" ]; then
	$FLOW  -i input -o output -s ./flow.con
	$FLOW  -i input -o output -s ./flow_vtk.con
	$FLOW  -i input -o output -s ./flow_mumps.con
	$FLOW  -i input -o output -s ./flow_matis.con
else
	$FLOW  -i input -o output -s ./$CONFIG.con	#-ksp_atol 1.0e-10 -ksp_rtol 1.0e-10 -ksp_monitor
fi


echo
echo "################################################################################"
echo "##########                     01_steady_flow_132d                   ##  End  ##"
echo "################################################################################"
echo


