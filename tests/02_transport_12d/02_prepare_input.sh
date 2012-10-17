#! /bin/bash	
echo
echo "##########               02_transport_12d: NGH and BCD               ## Start ##"
echo

# paths to trunk/bin to NGH and BCD programs
NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

#NAME="test1"

#meshing geo file
#gmsh -2 -algo front2d -3 -algo front3d -o input/$NAME.msh -format msh $NAME.geo

# generating neighbouring
$NGH ./ngh.ini

# preparing boundary condition files (flow, transport)
$BCD ./bcd.ini

echo
echo "##########               02_transport_12d: NGH and BCD               ### End ###"
echo


