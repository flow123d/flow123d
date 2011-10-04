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
#set -x
# This script compares two given diretories

#set if nor running as part of run_test.sh script
#TEST_DIR=/cygdrive/c/cygwin_drive/flow/branches/1.6.5/tests/

# Get actual output dir as the first argument
OUT="$1"

# directory with reference output
REF_OUT="$2"

# optional arguments to proper displaying
# ini file
INI=$3
# num. of  procs
n=$4

# file comparison script, if not running as part of run_test.sh, set absolute path to ndiff
# ex. NDIFF="/cygdrive/c/cygwin_drive/flow/branches/1.6.5/tests/ndiff.pl"
NDIFF="${TEST_DIR}/../ndiff.pl"

if [ ! -x "$NDIFF" ]
then 
  echo "can not find or run ndiff.pl"
  exit 1
fi

INI="${INI##*/}"
TEST="${TEST_DIR##*/}"
SCRIPT_PATH_DIR="`pwd`"
ERROR=0


if ! cd "$REF_OUT" 
then
  echo "Error: can not change to directory with reference output!"
  exit 1 
fi
	
#VAR which contains names of files
VAR=`ls`
	
#cycle for checking output files
for x in $VAR
do
if [ -d "$x" ]; then
	if ! "$SCRIPT_PATH_DIR/../run_check.sh" "$OUT/$x/" "$REF_OUT/$x/"; then
		ERROR=1
		echo "ERROR"
	fi	
else
	if [ $x == err ]; then
		echo "" | tee --append "${TEST_DIR}/stdout_diff.log"
		echo "Err log : ini file: ${INI}, procs: ${n}, test: ${TEST}" | tee --append "${TEST_DIR}/stdout_diff.log"
		touch empty
		if ! "$NDIFF" -o "${TEST_DIR}/diff.log" "$OUT/$x" empty | tee --append "${TEST_DIR}/stdout_diff.log"; then
			ERROR=1
		fi
		rm empty
	#else compare rest of files
	elif [ $x == out ]; then
		echo ""
	else 
		if [ -a "$OUT/$x" ]; then			
			echo "" | tee --append "${TEST_DIR}/stdout_diff.log"
			echo "Test: ${TEST}, ini file: ${INI}, procs: ${n}, file: $x " | tee --append "${TEST_DIR}/stdout_diff.log"
			if ! "$NDIFF" -o "${TEST_DIR}/diff.log" "$OUT/$x" "$REF_OUT/$x" | tee --append "${TEST_DIR}/stdout_diff.log"; then
				ERROR=1
			fi
		else
			echo "Error: Missing output file: $x" | tee --append "${TEST_DIR}/stdout_diff.log"
			ERROR=1
		fi
	fi
fi
done

if [ $ERROR == 1 ]; then
	exit 1
fi

