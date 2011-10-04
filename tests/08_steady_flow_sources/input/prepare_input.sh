set -x

FLOW=$HOME/workspace/flow123d
NGH=$FLOW/external_software/ngh/bin/ngh
BCD=$FLOW/external_software/bcd/bin/bcd

# assume uloha.geo
gmsh -2 test008.geo
# generate neighbouring
$NGH ngh.ini

#prepare  transport variable boundary 
$BCD "bcd.ini"

