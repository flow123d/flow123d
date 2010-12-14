#!/bin/bash

# passing arguments
while [ \( -n "$1" \) -a \( ! "$1" == "--" \) ]
do
  if [ "$1" == "-np" ]; then
    shift
	NPROC=$1				
	shift
  elif [ "$1" == "-ini" ]; then
    shift
    INI_FILE=$1
	export INI=${INI_FILE##*/}
	SOURCE_DIR=${INI_FILE%/*}
  elif [ "$1" == "--flow_params" ]; then
    shift
    FLOW_PARAMS=$1
  elif [ "$1" == "-h" ]; then
    echo " This is Flow123d help page:
	args:
	-np 		set number of procs
	-ini 		set absolut or relative path to ini file
	-wtime 		set maximal time to wait to finish job"
	break
    shift
  else
    shift
  fi
done 

# set path to script dir
export SCRIPT_PATH_DIR=$PWD

# set path to test dir
#FILE_PATH_DIR=$SOURCE_DIR

# path to MPIEXEC
MPI_RUN=$SCRIPT_PATH_DIR/../bin/mpiexec

#paths to dirs, relative to tests/$WRK
if [ -e $SCRIPT_PATH_DIR/../bin/flow123d.exe ]; then
	export EXECUTABLE=$SCRIPT_PATH_DIR/../bin/flow123d.exe;
elif [ -e $SCRIPT_PATH_DIR/../bin/flow123d ]; then
	export EXECUTABLE=$SCRIPT_PATH_DIR/../bin/flow123d;
else
	echo "Error: missing executable file!"
	exit 1
fi
 
 cd $SOURCE_DIR

#run flow, check if exists mpiexec skript, else allow run only with 1 procs without MPIEXEC	
if [ -n "$MACHINE_NAME" ]; then
	../${MACHINE_NAME}_make_pbs.sh
	qsub -pe mpi $NPROC ${MACHINE_NAME}_run_pbs.qsub
else
	if [ -e $MPI_RUN ]; then
		$MPI_RUN -np $NPROC $EXECUTABLE -S $INI $FLOW_PARAMS 2>err 1>out
	else 
		echo "Error: Missing mpiexec, unavailable to proceed with more then one procs"
		exit 1
	fi
fi
	
	