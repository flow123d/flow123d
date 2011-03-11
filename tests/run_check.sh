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

# Get actual output dir as the first argument
OUT=$1

REF_OUT=$2

# file comparison script
NDIFF="`pwd`/../ndiff.pl"
if [ ! -x "$NDIFF" ]
then echo "can not find or run ndiff.pl"
fi

#names of ini files
INI_FILE=$1

INI=${INI_FILE##*/}
SOURCE_DIR=${INI_FILE%/*}

#number of processors

SCRIPT_PATH_DIR="$PWD"

ERROR=0

suff=.ini
	
cd "$REF_OUT"
	
#VAR which contains names of files
VAR=`ls|sort`
	
#cycle for checking output files
for x in $VAR
do
#if its err file, compare with empty file
if [ "$x" == "err" ]; then
	echo "" | tee stdout_diff.log
	echo "Err log:" | tee stdout_diff.log
	touch empty
	if ! "$NDIFF" -o diff.log "$OUT/$x" empty; then
		ERROR=1
	fi
	rm empty
#else compare rest of files
elif [ "$x" == "out" ]; then
	echo ""
else 
	if [ -a "$OUT/$x" ]; then			
		echo "" | tee stdout_diff.log
		echo "File: $x" | tee -a stdout_diff.log
		if ! "$NDIFF" -o diff.log "$OUT/$x" "$REF_OUT/$x"; then
			ERROR=1
		fi
	else
		echo "Error: Missing output file: $x" | tee -a stdout_diff.log
		ERROR=1
	fi
fi
done

mv diff.log "$OUT"
mv stdout_diff.log "$OUT"

if [ $ERROR == 1 ]; then
	exit 1
fi

