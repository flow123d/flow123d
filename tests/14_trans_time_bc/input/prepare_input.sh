set -x

FLOW=$HOME/workspace/flow123d
NGH=$FLOW/external_software/ngh/bin/ngh
BCD=$FLOW/external_software/bcd/bin/bcd

# assume uloha.geo
gmsh -2 uloha.geo
# generate neighbouring
$NGH ngh.ini

#prepare  transport variable boundary 
$BCD "bcd.ini"

tlevel=0
for value in 20 0 40 0
do
  new_name="uloha_tbc_`printf %03d $tlevel`"
  sed "s/0\\.000000000000/$value/" <uloha_tbc >$new_name
  ((tlevel++))
done