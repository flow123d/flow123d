#!/bin/bash

#names of ini files
export INI_FILE=$1

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

# set path to script dir
#export SCRIPT_PATH_DIR=$PWD

# set path to test dir
#export FILE_PATH_DIR=./../tests/01_steady_flow_123d/

# path to MPIEXEC
export MPI_RUN=$SCRIPT_PATH_DIR/../bin/mpiexec

#paths to dirs, relative to tests/$WRK
if [ -e $SCRIPT_PATH_DIR/../bin/flow123d.exe ]; then
	export EXECUTABLE=$SCRIPT_PATH_DIR/../bin/flow123d.exe;
elif [ -e $SCRIPT_PATH_DIR/../bin/flow123d ]; then
	export EXECUTABLE=$SCRIPT_PATH_DIR/../bin/flow123d;
else
	echo "Error: missing executable file!"
	exit 1
fi
 
 #SUFF=.ini
 #INI=$1
 
 cd $FILE_PATH_DIR
	
	#mkdir -p ./$INPUT
	#mkdir ./$WRK/$OUTPUT

	#copy input and ini file
	#cp ./$3/$INPUT/* ./$WRK/$INPUT
	#cp ./$3/$1 ./$WRK/

	#run flow, check if exists mpiexec skript, else allow run only with 1 procs without MPIEXEC
	#cd $WRK
	if [ -n "$MACHINE_NAME" ]; then
		../${MACHINE_NAME}_make_pbs.sh
		qsub -pe mpi $NPROC ${MACHINE_NAME}_run_pbs.qsub
	else
		if [ $2 -eq 1 ]; then
			$EXECUTABLE -S $INI_FILE $FLOW_PARAMS 2>err 1>out
		elif [ -e $MPI_RUN ]; then
			$MPI_RUN -np $2 $EXECUTABLE -S $INI_FILE $FLOW_PARAMS 2>err 1>out
		else 
			echo "Error: Missing mpiexec, unavailable to proceed with more then one procs"
			#cd ..
			#rm -rf $WRK
			exit 1
		fi
	fi
	
	#makes dir for results, copy output files and logs
	#if [ -d ../$RES/${INI%$SUFF}.$2/ ]; then
	#	rm -rf ../$RES/${INI%$SUFF}.$2/
	#	mkdir -p ../$RES/${INI%$SUFF}.$2/
	#else 
	#	mkdir -p ../$RES/${INI%$SUFF}.$2/
	#fi
	
	#cp ./err ../$RES/${INI%$SUFF}.$2/
	#cp ./out ../$RES/${INI%$SUFF}.$2/
	
	#VAR=`ls $FILE_PATH_DIR/output/`
	
	#if [ -z "$VAR" ]; then
	#	echo "Error: Missing output files"
		#cd ..
		#rm -rf $WRK
	#	exit 1
	#else
	#	cp ./output/* ../$RES/${INI%$SUFF}.$2/
	#fi
	
	#remove tmp dir
	#cd ..
	#rm -rf $WRK