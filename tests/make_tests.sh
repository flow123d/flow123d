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

#names of tests
SOURCE_DIR=$3

#additional flow parameters
FLOW_PARAMS=$4

#input dir
INPUT=input

#result dir
RES=Results

#input dir
OUTPUT=output

#paths to dirs, relative to tests/$WRK
if [ -e ../bin/flow123d.exe ]; then
	EXECUTABLE=../../bin/flow123d.exe;
elif [ -e ../bin/flow123d ]; then
	EXECUTABLE=../../bin/flow123d;
else
	echo "Error: missing executable file!"
	exit 1
fi

WRK=tmp
MPI_RUN=../../bin/mpiexec

#tests of scripts permission
if [ -x ndiff.pl ]; then
	echo "ndiff.pl: permission ok"
else 
	chmod u+x ndiff.pl
fi 

if [ -x make_tests.sh ]; then	
	echo "make_tests.sh: permission ok"
else 
	chmod u+x make_tests.sh
 fi
 
 #script for executing flow123d and copying results
function run_flow() {
	#$1 - name of ini file
	#$2 - nprocs
	#$3 - source_dir
		
	SUFF=.ini
	INI=$1
	
	mkdir -p ./$WRK/$INPUT
	mkdir ./$WRK/$OUTPUT

	#copy input and ini file
	cp ./$3/$INPUT/* ./$WRK/$INPUT
	cp ./$3/$1 ./$WRK/

	#run flow, check if exists mpiexec skript, else allow run only with 1 procs without MPIEXEC
	cd $WRK
	if [ $2 -eq 1 ]; then
		$EXECUTABLE -s $1 $FLOW_PARAMS 2>err 1>out
	elif [ -e $MPI_RUN ]; then
		$MPI_RUN -np $2 $EXECUTABLE -s $1 $FLOW_PARAMS 2>err 1>out
	else 
		echo "Error: Missing mpiexec, unavailable to proceed with more then one procs"
		cd ..
		rm -rf $WRK
		exit 1
	fi
	
	#makes dir for results, copy output files and logs
	if [ -d ../$3/$RES/${INI%$SUFF}.$2/ ]; then
		rm -rf ../$3/$RES/${INI%$SUFF}.$2/
		mkdir -p ../$3/$RES/${INI%$SUFF}.$2/
	else 
		mkdir -p ../$3/$RES/${INI%$SUFF}.$2/
	fi
	
	cp ./err ../$3/$RES/${INI%$SUFF}.$2/
	cp ./out ../$3/$RES/${INI%$SUFF}.$2/
	
	VAR=`ls ./$OUTPUT/`
	
	if [ -z "$VAR" ]; then
		echo "Error: Missing output files"
		cd ..
		rm -rf $WRK
		exit 1
	else
		cp ./$OUTPUT/* ../$3/$RES/${INI%$SUFF}.$2/
	fi
	
	#remove tmp dir
	cd ..
	rm -rf $WRK
	
}

#script for checking new output and reference output
function check() {
	
	#suffix and name of inifile
	suff=.ini
	ini=$1
	
	cd $3/$RES/${ini%$suff}.$2
	
	#VAR which contains names of files
	VAR=`ls|sort`
	cd ../..
	
	#cycle for checking output files
	for x in $VAR
	do
	#if its err file, compare with empty file
	if [ "$x" == "err" ]; then
		cd ..
		echo "" | tee stdout_diff.log
		echo "Err log:" | tee stdout_diff.log
		touch empty
		./ndiff.pl -o diff.log ./$3/$RES/${ini%$suff}.$2/$x empty | tee -a stdout_diff.log
		rm empty
		cd $3
	#else compare rest of files
	elif [ "$x" == "out" ]; then
		touch empty
		rm empty
	else 
		cd ..
		if [ -a ./$3/$RES/${ini%$suff}.$2/$x ]; then
			
			echo "" | tee stdout_diff.log
			echo "Test:$3, ini file: $1, n proc: $2, output file: $x" | tee -a stdout_diff.log
			./ndiff.pl -o diff.log ./$3/$RES/${ini%$suff}.$2/$x ./$3/output/$x | tee -a stdout_diff.log
			cd $3
		else
			echo "Error: Missing one or more output file.($x)"
			exit 1
		fi
	fi
	done
	mv ../diff.log ./$RES/${ini%$suff}.$2/
	mv ../stdout_diff.log ./$RES/${ini%$suff}.$2/
	cd ..
}

#cycles for more ini files and nproc
for f in $INI_FILE
do
for n in $NPROC
do
for t in $SOURCE_DIR
do
#runs script which copies flow input files and ini files from given folders
run_flow $f $n $t

#runs ndiff.pl skript with ref and computed output files
echo "******************************************"
check $f $n $t
echo "******************************************"

done
done
done