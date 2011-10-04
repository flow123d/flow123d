set -x

FLOW=$HOME/workspace/flow123d
NGH=$FLOW/external_software/ngh/bin/ngh
BCD=$FLOW/external_software/bcd/bin/bcd

# assume uloha.geo
# algo: meshadapt, del2d, front2d
gmsh -2 -algo del2d uloha.geo
# generate neighbouring
$NGH ngh.ini

#prepare  transport variable boundary 
$BCD "bcd.ini"

# use tic for initial pressure
cat uloha.tic | sed 's/Concentrations/PressureHead/' | sed 's/^\([0-9]*\)[ \t]*[0-9]*/\1/' > pressure_initial.in


tlevel=0
for value in 10 0
do
  new_name="uloha_tbc_`printf %03d $tlevel`"
  sed "s/10\\.0*/$value/" <uloha_tbc >$new_name
  ((tlevel++))
done