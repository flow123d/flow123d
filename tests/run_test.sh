#!/bin/bash

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
$SCRIPT_PATH_DIR/run_flow.sh -ini $INI_FILE -np $n --flow_params $FLOW_PARAMS

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