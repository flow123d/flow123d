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

#names of ini files
INI_FILE=$1

#number of processors
NPROC=$2

#adition flow params
FLOW_PARAMS=$3

# set path to script dir
SCRIPT_PATH_DIR=$PWD

# set path to test dir
#PATH_DIR=/c/cyg_drive/tests/bla.txt

INI=${INI_FILE##*/}
SOURCE_DIR=${INI_FILE%/*}

SUFF=.ini

for n in $NPROC
do
#runs script which copies flow input files and ini files from given folders
export FILE_PATH_DIR=$SOURCE_DIR
$SCRIPT_PATH_DIR/run_flow.sh -ini $INI_FILE -np $n -- $FLOW_PARAMS

#while [ ! -f $FILE_PATH_DIR/out ]
#do	
#	while [ -e $FILE_PATH_DIR/lock ]
#	do
#		sleep 20
#	done
#done

if [ -d $FILE_PATH_DIR/Results/${INI%$SUFF}.$2/ ]; then
		rm -rf $FILE_PATH_DIR/Results/${INI%$SUFF}.$2/
		mkdir -p $FILE_PATH_DIR/Results/${INI%$SUFF}.$2/
	else 
		mkdir -p $FILE_PATH_DIR/Results/${INI%$SUFF}.$2/
	fi

mv $FILE_PATH_DIR/err $FILE_PATH_DIR/Results/${INI%$SUFF}.$2/
mv $FILE_PATH_DIR/out $FILE_PATH_DIR/Results/${INI%$SUFF}.$2/

#runs ndiff.pl skript with ref and computed output files
echo "******************************************"
$SCRIPT_PATH_DIR/run_check.sh $INI_FILE $n
echo "******************************************"

done