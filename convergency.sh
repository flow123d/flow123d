#!/bin/bash

# script for running solution of richards problem for various values of discretization parameters
# in order to produce convergency graph

#set -x

MAIN_DIR=conv_results
mkdir $MAIN_DIR

finest=$MAIN_DIR/dh_0.0001_dt_0.1
echo "time        dh         dt          p_err            q_err" >$MAIN_DIR/summary
for dh in 0.1 0.05 0.02 0.01 0.005 0.002 0.001 0.0005 0.0002 
do
  for dt in 0.1 50 20 10 5.0 2.0 1.0 0.5 0.2
 #for dt in 0.2
  do
    #dt=$dh
    res_dir="$MAIN_DIR/dh_${dh}_dt_${dt}"
    echo " computing $res_dir"  
    if [ -z "$finest" ]
    then
      finest=$res_dir
    fi
    echo "fin: $finest"  

    mkdir $res_dir
    rm -f output/*
    cat input.template.txt |sed "s/DH/$dh/" | sed "s/DT/$dt/" > input.txt
    /usr/bin/time -a -o out ./lib/richards-2d >out 
    grep "user" out 

    cp -r output $res_dir
    cp out $res_dir
    
    grep "PRINT" out > errors
    last=`tail -n 1 errors | sed "s/.*: //"`
    time=${last}
    #both_err=`echo $last | sed "s/.*: //"`
    #p_err=${both_err%% *}
    #q_err=${both_err##* }
    #mv errors $res_dir

    
    num_p_err=`~/local/ParaView-3.12.0-Linux-x86_64/bin/pvpython \
    pv_compute.py $res_dir/output/solution-020.vtk $finest/output/solution-020.vtk post_phead \
    | grep 'L2' | sed 's/^.*: //'` 

    num_q_err=`~/local/ParaView-3.12.0-Linux-x86_64/bin/pvpython \
    pv_compute.py $res_dir/output/solution-020.vtk $finest/output/solution-020.vtk flux \
    | grep 'L2' | sed 's/^.*: //' `
    
    echo "$time  $dh    $dt    $num_p_err $num_q_err"
    echo "$time  $dh    $dt    $num_p_err $num_q_err" >>$MAIN_DIR/summary
    

  done
done
    