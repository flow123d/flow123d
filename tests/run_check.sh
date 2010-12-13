#!/bin/bash

#names of ini files
INI_FILE=$1

#number of processors
NPROC=$2

#names of tests
SOURCE_DIR=$3

#input dir
INPUT=input

#result dir
RES=Results

#input dir
OUTPUT=output

if [ -x ndiff.pl ]; then
	echo "ndiff.pl: permission ok"
else 
	chmod u+x ndiff.pl
fi

suff=.ini
	ini=$1
	
	cd $FILE_PATH_DIR/ref_output
	
	#VAR which contains names of files
	VAR=`ls|sort`
	cd ..
	
	#cycle for checking output files
	for x in $VAR
	do
	#if its err file, compare with empty file
	if [ "$x" == "err" ]; then
		echo "" | tee stdout_diff.log
		echo "Err log:" | tee stdout_diff.log
		touch empty
		$SCRIPT_PATH_DIR/ndiff.pl -o diff.log $FILE_PATH_DIR/output/$x empty | tee -a stdout_diff.log
		rm empty
	#else compare rest of files
	elif [ "$x" == "out" ]; then
		echo ""
	else 
		if [ -a $FILE_PATH_DIR/output/$x ]; then			
			echo "" | tee stdout_diff.log
			echo "Test:$3, ini file: $1, n proc: $2, output file: $x" | tee -a stdout_diff.log
			$SCRIPT_PATH_DIR/ndiff.pl -o diff.log $FILE_PATH_DIR/output/$x $FILE_PATH_DIR/ref_output/$x | tee -a stdout_diff.log
		else
			echo "Error: Missing one or more output file.($x)"
			exit 1
		fi
	fi
	done
	mv $FILE_PATH_DIR/diff.log $FILE_PATH_DIR/Results/${ini%$suff}.$2/
	mv $FILE_PATH_DIR/stdout_diff.log $FILE_PATH_DIR/Results/${ini%$suff}.$2/