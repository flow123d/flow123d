#!/bin/bash

# script for running solution of richards problem for various values of discretization parameters
# in order to produce convergency graph

set -x

MAIN_DIR=conv_results
mkdir $MAIN_DIR

finest_solution=$MAIN_DIR/cn_trap_dh_0.001_dt_0.5

compare_times="50 100 200 400 800 1000"  
dh_steps="0.0002"
#dt_steps="0.5 0.2 0.1 0.05 0.02 0.01"
dt_steps="0.1"

compare_only=$1

echo "time        dh         dt          p_err            q_err         q_half_err" >$MAIN_DIR/summary

# path to Paraview python
pv_python="$HOME/local/ParaView-3.12.0-Linux-x86_64/bin/pvpython"

# path to python script to compare difference in of two VTK files
pv_compute="pv_compute.py"

# prepare input, compute, move results
# 
# * make input.txt from input.template.txt
# * run simulation
# * save results into MAIN_DIR
function compute {
    dh=$1
    dt=$2
    res_dir=$3

    echo " computing $res_dir"  

    mkdir $res_dir
    rm -f output/*
    cat input.template.txt |sed "s/DH/$dh/" | sed "s/DT/$dt/" > input.txt
    /usr/bin/time -a -o out ./lib/richards-2d >out 
    grep "user" out 

    mkdir $res_dir/output
    mv output/* $res_dir/output
    cat out |grep PRINT >$res_dir/print_times
    mv out $res_dir
    mv bc_output.out $res_dir
}

# extract errors computed by simulator itself and output them into summary file
# NOT TESTED YET
function compare_anal {
    dh=$1
    dt=$2
    res_dir=$3
    summary_file=$4

    # extract errors from stdout
    grep "L2" $res_dir/out > $res_dir/errors
    values=(`tail -n 1 $res_dir/errors | sed 's/Time: \([^ ]*\) L2 error: \(.*\)$/\1 \2/' `) # make array of 3 values: time, p_err, q_err
    time=${values[0]}
    p_err=${values[1]}
    q_err=${values[2]}

    echo "$time  $dh    $dt    $p_err $q_err"
    echo "$time  $dh    $dt    $p_err $q_err" >>$summary

}

# find output frame for given time and print_times file
function find_print_time {
  time=$1
  file=$2
  if [ ! -f $file ];then return; fi 
  cat $file | perl -e \
    "\$err=1000;
     while (\$_=<>) {s/PRINT time \(([0-9]*)\): *([^ ]*).*/\2 \1/; 
       if (abs( \$2 - $time) < \$err) 
       { \$err=abs(\$2 - $time); \$best_t=\$2; \$best_frame=\$1; 
     #print \"\$err \$best_t \$best_frame\\n\";
     }};
     if (\$err < \$best_t * 0.001) { print \$best_frame; }"
}

function calc {
bc << EOF
scale=8
$@
quit
EOF
}

function compare_num {
    local dh=$1
    local dt=$2
    local res_dir=$3
    local summary_file=$4

    # use variable compare_times and finest_solution
    local finest_dh=`get_finest_dh $finest_solution`  
    local finest_dt=`get_finest_dt $finest_solution`

    # compute error numericaly
    for time in $compare_times
    do
      half_time=`calc "$time - $dt / 2"`
      fin_sol_frame=`find_print_time $time $finest_solution/print_times`
      fin_sol_frame=`printf "solution-%03d.vtk" $fin_sol_frame`

      sol_frame=`find_print_time $time $res_dir/print_times`
      if [ -z $sol_frame ]
      then
        continue
      fi
      sol_frame=`printf "solution-%03d.vtk" $sol_frame`

      fin_sol_half_frame=`find_print_time $half_time $finest_solution/print_times`
      fin_sol_half_frame=`printf "solution-%03d.vtk" $fin_sol_half_frame`

      p_err=`$pv_python $pv_compute $res_dir/output/$sol_frame $finest_solution/output/$fin_sol_frame post_phead | grep 'L2' | sed 's/^.*: //'` 

      q_err=`$pv_python $pv_compute $res_dir/output/$sol_frame $finest_solution/output/$fin_sol_frame flux | grep 'L2' | sed 's/^.*: //'` 

      q_half_err=`$pv_python $pv_compute $res_dir/output/$sol_frame $finest_solution/output/$fin_sol_half_frame half_flux | grep 'L2' | sed 's/^.*: //'` 

    echo "$time  $dh    $dt    $p_err $q_err $q_half_err"
    echo "$time  $dh    $dt    $p_err $q_err $q_half_err" >>$summary_file

    done

}

function get_finest_dh {
  echo $1 | sed 's/.*dh_\([^_]*\).*/\1/'
}

function get_finest_dt {
  echo $1 | sed 's/.*dt_\([^_]*\).*/\1/'
}


function compute_finest {
  fin_sol=$1
}


###################################################33 MAIN CYCLE
#find_print_time 500.3 conv_results/dh_0.001_dt_100/print_times
#exit   

for dh in $dh_steps
do
  for dt in $dt_steps
  do
    res_dir="$MAIN_DIR/dh_${dh}_dt_${dt}"

    if [ ! "$compare_only" == "co" ] 
    then
      compute $dh $dt $res_dir
    fi
    
    # compare_anal $dh $dt $res_dir "$MAIN_DIR/summary"

    compare_num $dh $dt $res_dir "$MAIN_DIR/summary"  

  done
done
    