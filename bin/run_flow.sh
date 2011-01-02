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

set -x 

NPROC=1
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
  elif [ "$1" == "-q" ]; then
    shift
    QueueTime=$1
  elif [ "$1" == "-h" ]; then
    echo " This is Flow123d help page:
	args:
	-np 		set number of procs
	-ini 		set absolut or relative path to ini file
	-q 		set maximal time to wait to finish job"
	break
    shift
  else
    shift
  fi
done 

FLOW_PARAMS="$@"

if [ -z $INI_FILE ] 
then
  echo "Error ..."
  exit 1
fi

export INI=${INI_FILE##*/}
if [ "${INI_FILE%%[^/]*}" == "" ]
then 
  # relative path
  INI_FILE="/$INI_FILE"
  echo ${INI_FILE%/*}
  SOURCE_DIR="`pwd`/${INI_FILE%/*}"
else
  #absolute path
  SOURCE_DIR=${INI_FILE%/*}
fi


# set path to script dir + exports for make_pbs scripts
# TODO: co kdyz bude skript volan s absolutini cestou !!
export SCRIPT_PATH_DIR="`pwd`/${0%/*}" 

# path to MPIEXEC
MPI_RUN=$SCRIPT_PATH_DIR/mpiexec

#paths to dirs, relative to tests/$WRK
if [ -e $SCRIPT_PATH_DIR/flow123d.exe ]; then
	EXECUTABLE=$SCRIPT_PATH_DIR/flow123d.exe;
elif [ -e $SCRIPT_PATH_DIR/flow123d ]; then
	EXECUTABLE=$SCRIPT_PATH_DIR/flow123d;
else
	echo "Error: missing executable file!"
	exit 1
fi
 
#exports for make_pbs scripts
export MPI_RUN
export EXECUTABLE


cd $SOURCE_DIR

#run flow, check if exists mpiexec skript, else allow run only with 1 procs without MPIEXEC	
if [ -n "$MACHINE_NAME" ]; then
	../${MACHINE_NAME}_make_pbs.sh
	qsub -pe orte $NPROC ${MACHINE_NAME}_run_pbs.qsub
else
	if [ -e $MPI_RUN ]; then
		$MPI_RUN -np $NPROC $EXECUTABLE -S $INI $FLOW_PARAMS 2>err 1>out
	else 
		echo "Error: Missing mpiexec, unavailable to proceed with more then one procs"
		exit 1
	fi
fi
	
	