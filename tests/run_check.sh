#!/bin/bash

#names of ini files
INI_FILE=$1

INI=${INI_FILE##*/}
SOURCE_DIR=${INI_FILE%/*}

#number of processors
NPROC=$2

SCRIPT_PATH_DIR=$PWD

if [ -x ndiff.pl ]; then
	echo "ndiff.pl: permission ok"
else 
	chmod u+x ndiff.pl
fi

suff=.ini
	
	cd $SOURCE_DIR/ref_output
	
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
		$SCRIPT_PATH_DIR/ndiff.pl -o diff.log $SOURCE_DIR/output/$x empty | tee -a stdout_diff.log
		rm empty
	#else compare rest of files
	elif [ "$x" == "out" ]; then
		echo ""
	else 
		if [ -a $SOURCE_DIR/output/$x ]; then			
			echo "" | tee stdout_diff.log
			echo "Test:$3, ini file: $1, n proc: $2, output file: $x" | tee -a stdout_diff.log
			$SCRIPT_PATH_DIR/ndiff.pl -o diff.log $SOURCE_DIR/output/$x $SOURCE_DIR/ref_output/$x | tee -a stdout_diff.log
		else
			echo "Error: Missing one or more output file.($x)"
			exit 1
		fi
	fi
	done
	mv $SOURCE_DIR/diff.log $SOURCE_DIR/Results/${INI%$suff}.$2/
	mv $SOURCE_DIR/stdout_diff.log $SOURCE_DIR/Results/${INI%$suff}.$2/