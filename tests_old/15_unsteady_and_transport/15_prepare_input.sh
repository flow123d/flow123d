#! /bin/bash	
echo
echo "##########        15_unsteady_transport: MESH and NGH and BCD        ### Start ###"
echo

NGH=../../bin/ngh/bin/ngh
BCD=../../bin/bcd/bin/bcd

NAME="test15"

#generates a little bit different mesh then the original one, 
#meshing geo file
# algo: meshadapt, del2d, front2d
#gmsh -2 -algo del2d -o input/$NAME.msh $NAME.geo


$NGH ./ngh.ini
$BCD ./bcd.ini


# use tic file to make initial pressure file
cat input/test15.tic | sed 's/Concentrations/PressureHead/' > input/pressure_initial.in

#making time variable bc for transport from the one file (replacing values)
tlevel=0
for value in 10 0
do
  new_name="input/test15_tbc_`printf %03d $tlevel`"
  sed "s/0\\.000000000000/$value/" <input/test15_tbc >$new_name
  ((tlevel++))
done

#possible deleting original tbc file
rm input/test15_tbc

echo
echo "##########        15_unsteady_transport: MESH and NGH and BCD          ### End ###"
echo


