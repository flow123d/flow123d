#! /bin/bash	
echo
echo "##########         14_trans_time_bc: MESH and NGH and BCD          ### Start ###"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test14"

#generates a little bit different mesh then the original one, 
#meshing geo file
#gmsh -2 -o input/$NAME.msh -format msh $NAME.geo

$NGH ./ngh.ini
$BCD ./bcd.ini


#making time variable bc for transport from the one file (replacing values)
tlevel=0
for value in 20 0 40 0
do
  new_name="input/test14_tbc_`printf %03d $tlevel`"
  sed "s/0\\.000000000000/$value/" <input/test14.tbc >$new_name
  ((tlevel++))
done

rm input/test14.tbc

echo
echo "##########         14_trans_time_bc: MESH and NGH and BCD          ### End ###"
echo


