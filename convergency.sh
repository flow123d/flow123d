#!/bin/bash

# script for running solution of richards problem for various values of discretization parameters
# in order to produce convergency graph

#set -x

MAIN_DIR=conv_thin
mkdir $MAIN_DIR

echo "time        dh         dt          p_err            q_err" >$MAIN_DIR/summary
for dh in 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001
do
  for dt in 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001
  do
    res_dir="$MAIN_DIR/dh_${dh}_dt_${dt}"
    echo " computing $res_dir"  
    mkdir $res_dir
    cat input.template.txt |sed "s/DH/$dh/" | sed "s/DT/$dt/" > input.txt
    ./lib/richards-2d > out
    grep "L2" out > errors
    last=`tail -n 1 errors | sed "s/Time: //"`
    time=${last%%L2*}
    both_err=`echo $last | sed "s/.*: //"`
    p_err=${both_err%% *}
    q_err=${both_err##* }
    entry="$time  $dh    $dt     $p_err         $q_err"
    echo $entry >>$MAIN_DIR/summary
    mv out $res_dir
    mv errors $res_dir
  done
done
    