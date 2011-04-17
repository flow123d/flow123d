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
set -x
# This script assumes that it is running in particular subdir of the "tests"
# dir. 
pwd
export TEST_DIR="`pwd`"

#name of ini file
INI_FILES="$1"

#numbers of processors to run on
export NPROC="$2"

#adition flow params
FLOW_PARAMS="$3"

# how to run flow
RUN_FLOW=../../bin/run_flow.sh

ERROR=0

for INI_FILE in $INI_FILES
do
	for n in $NPROC
	do
		if ! $RUN_FLOW -s "$INI_FILE" -np "$n" -- "$FLOW_PARAMS"; then
			if [ ! -e ./out ]; then
				ERROR=1
				break 2
			else
				ERROR=1
				break
			fi
		fi

		if [ -e ./lock ]; then
			for i in $(seq 1 10)
			do
				if [! -e ./out ]; then
					sleep 30
				else
					break
				fi
			done
		
			if [! -e ./out ]; then
				echo "ERROR: Directory locked, no output file created, aborting"
				exit 1
			fi
		fi

		for i in $(seq 1 40)
		do	
			if [ -e ./lock ]; then
				sleep 30
			else 
				break
			fi
			if [ $i == 40 ]; then
				echo "Error, time run out, exit 1"
				exit 1
			fi
		done
			
		SAVE_OUTPUT="$TEST_DIR/Results/${INI_FILE%.ini}.$n"
		if [ -d "$SAVE_OUTPUT" ]; then
			rm -rf "$SAVE_OUTPUT"
			mkdir -p "$SAVE_OUTPUT"
		else 
			mkdir -p "$SAVE_OUTPUT"
		fi

		mv ./err "$SAVE_OUTPUT"
		mv ./out "$SAVE_OUTPUT"
		mv ./*.log "$SAVE_OUTPUT"
		mv ./output/* "$SAVE_OUTPUT"
		
		
		#runs ndiff.pl skript with ref and computed output files
		echo "******************************************"
		if ! ../run_check.sh "$SAVE_OUTPUT" "$TEST_DIR/ref_output" "$INI_FILE" "$n"; then
			ERROR=1
		fi
		echo "******************************************"
	
		mv "${TEST_DIR}/diff.log" "$SAVE_OUTPUT"
		mv "${TEST_DIR}/stdout_diff.log" "$SAVE_OUTPUT"
  
	done
done

if [ $ERROR == 1 ]; then
	exit 1
fi
