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

#check if script is called with relative or absolute path

ShowHelp() {
echo " This is Flow123d help page. Syntax:
 
  run_flow.sh -np N -s INI_FILE [-q TIME] [-m MACHINE]

	args:
	-np N		set number of procs N
	-s INI_FILE	set absolut or relative path to ini file
	-q TIME		set maximal TIME to wait to finish job
	-m MACHINE	name of the machine, to determine particular start script for PBS, 
			default can be specified in makefile.in"
}

if [ ! "${0%%[^/]*}" == "" ]; then
	SCRIPT_PATH_DIR="${0%/*}"
else	
	SCRIPT_PATH_DIR="`pwd`/${0%/*}"
fi	
MACHINE_SCRIPT="$SCRIPT_PATH_DIR/current_flow"
NPROC=1
# passing arguments
while [ \( -n "$1" \) -a \( ! "$1" == "--" \) ]
do
  if [ "$1" == "-np" ]; then
    shift
	NPROC=$1				
	shift
  elif [ "$1" == "-s" ]; then
    shift
    INI_FILE=$1
  elif [ "$1" == "-q" ]; then
    shift
    QueueTime=$1
  elif [ "$1" == "-m" ]; then
    shift
	if [ -e "$SCRIPT_PATH_DIR/${1}_flow.sh" ]; then
		MACHINE_SCRIPT="$SCRIPT_PATH_DIR/${1}_flow.sh"
	else
		echo "ERROR: Missing MACHINE script, using default."
	fi
  elif [ "$1" == "-h" ]; then
    ShowHelp
	break
    shift
  else
    shift
  fi
done 

FLOW_PARAMS="$@"

if [ -z "$INI_FILE" ] 
then
  echo "Error: No ini file!"
  ShowHelp
  exit 1
fi


export INI=${INI_FILE##*/}
if [ "${INI_FILE%%[^/]*}" == "" ]
then 
  # relative path
  INI_FILE="/$INI_FILE"
  SOURCE_DIR="`pwd`/${INI_FILE%/*}"
else
  #absolute path
  SOURCE_DIR="${INI_FILE%/*}"
fi



# path to MPIEXEC
MPI_RUN="$SCRIPT_PATH_DIR/mpiexec"

#paths to dirs, relative to tests/$WRK
if [ -e "$SCRIPT_PATH_DIR/flow123d.exe" ]; then
	EXECUTABLE="$SCRIPT_PATH_DIR/flow123d.exe";
elif [ -e "$SCRIPT_PATH_DIR/flow123d" ]; then
	EXECUTABLE="$SCRIPT_PATH_DIR/flow123d";
else
	echo "Error: missing executable file!"
	exit 1
fi
 
#exports for make_pbs scripts
export MPI_RUN
export EXECUTABLE
export NPROC
export INI
export FLOW_PARAMS
export SOURCE_DIR

cd "$SOURCE_DIR"

#run flow, check if exists mpiexec skript, else allow run only with 1 procs without MPIEXEC	

if [ -e "$MPI_RUN" ]; then
	#$MPI_RUN -np $NPROC $EXECUTABLE -s $INI $FLOW_PARAMS 2>err 1>out
	"$MACHINE_SCRIPT"
else 
	echo "Error: Missing mpiexec, unavailable to proceed with more then one procs"
	exit 1
fi


if [ -e ./lock ]; then
	for i in $(seq 1 10)
	do
		if [! -e ./out ]; then
			sleep 10
		else
			break
		fi
	done
	if [! -e ./out ]; then
		echo "ERROR: Directory locked, no output file created, aborting"
		exit 1
	fi
fi

for i in $(seq 1 10)
do	
	if [ -e ./lock ]; then
		sleep 10
	else 
		break
	fi
	if [ $i == 10 ]; then
		echo "Error, directory locked too long, exit 1"
		exit 1
	fi
done
	
	