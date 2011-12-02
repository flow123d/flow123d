#!/bin/bash

# script for printing 'summary' files produced by script convergency.sh
# 'summary' has on each line dx, dt, pressure_error, flux_error
#
# this script produce four graphs for dependency of pressure/flux error on dx/dt
# in each graph there are several lines for the second discretization parameter
#set -x

in_file=summary

# get possible dt and dx values
dt_values="`cat summary | grep -v "time" | sed "s/^2 //"| sed "s/[^ ]* //"|sed "s/ .*//"|sort -u `"
dx_values="`cat summary | grep -v "time" | sed "s/^2 //"| sed "s/ .*//"|sort -u `"


function plot_file {
  file_name=$1
  swap=$2
  out_var="$3"
  
  if [ "$swap" == "swap" ]
  then
    x_list="$dt_values"
    t_list="$dx_values"
  else
    x_list="$dx_values"
    t_list="$dt_values"
  fi

  echo >$file_name
  # dt on x axes
  for x in $x_list
  do
    line="$x "
    for t in $t_list
    do
      if [ "$swap" == "swap" ]
      then
        dx=$t; dt=$x
      else
        dx=$x;dt=$t
      fi
      
      q_err=`cat summary| grep "^2 $dx $dt" | sed "s/^2 [^ ]* [^ ]* [^ ]* //"`
      p_err=`cat summary| grep "^2 $dx $dt" | sed "s/^2 [^ ]* [^ ]* //" | sed "s/ .*//"`
      
      line+=" ${!out_var}"
    done
    echo $line >>$file_name
  done
}

function make_plot_cmd {
  file_name=$1
  x_axes_list=$2

  out="set title \"$file_name\";plot "
  i=2;
  for dx in $x_axes_list
  do
    out+="'$file_name' using 1:$i title '$dx', "
    i=`expr $i + 1`
  done
  echo ${out} 0.001*x with lines, 0.01*x**2 with lines
}

plot_file "p_err_on_dt" "swap" "p_err"
plot_file "q_err_on_dt" "swap" "q_err"
plot_file "p_err_on_dx" "no_swap" "p_err"
plot_file "q_err_on_dx" "no_swap" "q_err"


p_dt_plot=`make_plot_cmd p_err_on_dt "$dt_values"`
q_dt_plot=`make_plot_cmd q_err_on_dt "$dt_values"`
p_dx_plot=`make_plot_cmd p_err_on_dx "$dx_values"`
q_dx_plot=`make_plot_cmd q_err_on_dx "$dx_values"`

gnuplot <<END
set xr [0.001:0.5]
set yr [0.0001:0.1]
set logscale x
set logscale y
set style data linespoints
unset key

set size 1,1
set origin 0,0
set multiplot

# plot the first graph so that it takes a quarter of the screen
set size 0.5,0.5
set origin 0,0.5
$p_dt_plot

# plot the second graph so that it takes a quarter of the screen
set size 0.5,0.5
set origin 0,0
$q_dt_plot

# plot the third graph so that it takes a quarter of the screen
set size 0.5,0.5
set origin 0.5,0.5
$p_dx_plot

# plot the fourth graph so that it takes a quarter of the screen
set size 0.5,0.5
set origin 0.5,0
$q_dx_plot

# On some terminals, nothing gets plotted until this command is issued
unset multiplot
  
END