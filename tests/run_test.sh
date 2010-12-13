#!/bin/bash

#names of ini files
INI_FILE=$1

#number of processors
NPROC=$2

#names of tests
SOURCE_DIR=$3

# set path to script dir
export SCRIPT_PATH_DIR=$PWD

# set path to test dir
export FILE_PATH_DIR=$SCRIPT_PATH_DIR/../tests/

SUFF=.ini

for f in $INI_FILE
do
for n in $NPROC
do
for t in $SOURCE_DIR
do
#runs script which copies flow input files and ini files from given folders
export FILE_PATH_DIR=$SCRIPT_PATH_DIR/../tests/$t
$SCRIPT_PATH_DIR/run_flow.sh $f $n $t

while [ ! -f $FILE_PATH_DIR/out ]
do	
	while [ -e $FILE_PATH_DIR/lock ]
	do
		sleep 20
	done
done

if [ -d $FILE_PATH_DIR/Results/${INI_FILE%$SUFF}.$2/ ]; then
		rm -rf $FILE_PATH_DIR/Results/${INI_FILE%$SUFF}.$2/
		mkdir -p $FILE_PATH_DIR/Results/${INI_FILE%$SUFF}.$2/
	else 
		mkdir -p $FILE_PATH_DIR/Results/${INI_FILE%$SUFF}.$2/
	fi

mv $FILE_PATH_DIR/err $FILE_PATH_DIR/Results/${INI_FILE%$SUFF}.$2/
mv $FILE_PATH_DIR/out $FILE_PATH_DIR/Results/${INI_FILE%$SUFF}.$2/

#runs ndiff.pl skript with ref and computed output files
echo "******************************************"
$SCRIPT_PATH_DIR/run_check.sh $f $n $t
echo "******************************************"

done
done
done