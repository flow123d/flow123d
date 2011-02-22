#!/bin/bash
# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License version 3 as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not,
# write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#

# This script assumes that it is running in particular subdir of the "tests"
# dir. 
TEST_DIR=`pwd`

#name of ini file
INI_FILE=$1

#numbers of processors to run on
NPROC=$2

#adition flow params
FLOW_PARAMS=$3

# how to run flow
#RUN_FLOW=../../bin/run_flow.sh
RUN_FLOW=../bin/run_flow.sh

ERROR=0

for i in $INI_FILE
do
for n in $NPROC
do
  if ! $RUN_FLOW -s $INI_FILE -np $n -- $FLOW_PARAMS; then
	echo " Error occured during computation, leaving."
	exit 1
  fi

  while [ ! -e ./out ]
  do		
  	while [ -e ./lock ]
  	do
  		sleep 20
  	done
  done

  SAVE_OUTPUT="$TEST_DIR/Results/${INI_FILE%.ini}.$n"
  if [ -d "$SAVE_OUTPUT" ]; then
		rm -rf "$SAVE_OUTPUT"
		mkdir -p "$SAVE_OUTPUT"
  else 
		mkdir -p "$SAVE_OUTPUT"
  fi

  mv ./err $SAVE_OUTPUT
  mv ./out $SAVE_OUTPUT
  mv ./*.log $SAVE_OUTPUT
  mv ./output/* $SAVE_OUTPUT

  #runs ndiff.pl skript with ref and computed output files
  echo "******************************************"
  if ! ./run_check.sh $SAVE_OUTPUT "$TEST_DIR/ref_output"; then
	ERROR=1
  fi
  echo "******************************************"

done
done

if [ $ERROR == 1 ]; then
	exit 1
fi