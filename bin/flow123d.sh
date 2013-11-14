#!/bin/bash
# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Center: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
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
# $Id: run_test.sh 1270 2011-08-09 11:37:24Z jiri.hnidek $
# $Revision: 1270 $
# $LastChangedBy: jiri.hnidek $
# $LastChangedDate: 2011-08-09 13:37:24 +0200 (Út, 09 srp 2011) $
#
# Author(s): Jiri Hnidek <jiri.hnidek@tul.cz>
#
# TODO:
# alow timeout parameter in different units
# alow mem parameter in different units

# Uncomment following line, when you want to debug bash script
# set -x 

# Relative path to mpiexec from the directory, where this script is placed
MPIEXEC="./mpiexec"
# Relative path to mpiexec binary from current/working directory
MPIEXEC="${0%/*}/${MPIEXEC}"

TIME_LIMIT_SH="${0%/*}/time_limit.sh"

# Relative path to Flow123d binary from the directory,
# where this script is placed
FLOW123D="./flow123d"
# Relative or absolute path to Flow123d binary from current/working directory
FLOW123D_REL="${0%/*}/${FLOW123D}"
# Absolute path
FLOW123D=`which "${FLOW123D_REL}"`
# Absolute path to working directory
WORKDIR="${PWD}"

AWK="awk"

# Print help to this script
function print_help {
	echo "SYNTAX: flow123d.sh [OPTIONS] -- [FLOW_PARAMS]"
	echo ""
	echo "Submit a flow123d job to a batch system using particular backend script deduced from HOSTNAME or from the --host parameter."
	echo "If the backend script is not found for the host, the flow123d is run interactively."
	echo "All parameters after "--" will be passed to flow123d. "
	echo "Unresolved parameters will be passed to qsub command if it is used."
	echo ""
	echo "OPTIONS:"
	echo "    -h, --help                    Print this help"
	echo "    --host HOSTNAME               Use given HOSTNAME for backend script resolution. Script 'config/\${HOSTNAME}.sh' must exist."
	echo "    -t, --walltime TIMEOUT        Specifies a maximum time period after which Flow123d will be killed.  Time TIMEOUT is expressed in seconds as an integer,\n or in the form: [[hours:]minutes:]seconds[.milliseconds]."
	echo "    -m, --mem MEM                 Flow123d can use only MEM magabytes per process."
	echo "    -n NICE, --nice               Run Flow123d with changed (lower) priority."
	echo "    -np N                         Run Flow123d using N parallel processes." 
	echo "    -ppn PPN                      Run PPN processes per node. NP should be divisible by PPN other wise it will be truncated."
	echo "    -q, --queue QUEUE             Name of queue to use for batch processing. For interactive runs this redirect stdout and stderr to the file with name in format QUEUE.DATE."
	echo ""
#	echo "RETURN VALUES:"
#	echo "    0               Flow123d process exited normaly"
#	echo "    10              Bad argument/option"
#	echo "    11              Flow123d can't be executed"
#	echo "    12              INI file doesn't exist"
#	echo "    13              Flow123d can't read INI file"
#	echo "    14              Flow123d process was killed (TIMEOUT exceeded)"
#	echo "    15              Flow123d process crashed"
#	echo ""
}

# Parse script arguments
function parse_arguments()
{
	# Make sure that INI_FILE is not set
	#unset INI_FILE

	# Default number of processes
	NP=1

        # Defaulf backend hostname
        BACKEND_HOST=${HOSTNAME}
        
        # Default value for nice
        NICE=0

	# print help when called with no parameters
	if [ -z "$1" ]; then print_help; exit 1; fi
	
	# Parse arguments (do not use bash builtin command getopts, it is too restrictive)
	while [ -n "$1" ]
	do
		if [ "$1" == "-h" -o "$1" == "--help" ];
		then
			print_help;
			exit 0;
                elif [ "$1" == "--host" ];
                then
                        shift; BACKEND_HOST="$1"
                elif [ "$1" == "-t" -o "$1" == "--walltime" ];
                then
                        shift; TIMEOUT="$1"
                elif [ "$1" == "-m" -o "$1" == "--mem" ];
                then
                        shift; MEM="$1"
		elif [ "$1" == "-n" -o "$1" == "--nice" ];
		then
                        shift; NICE="$1"
                elif [ "$1" == "-p" -o "$1" == "-np" ];
                then
                        shift; NP="$1"
                elif [ "$1" == "-q" -o "$1" == "--queue" ];
                then
                        shift; QUEUE="$1"                     
		elif [ "$1" == "-ppn" ];
                then
                        shift; PPN="$1"
		elif [ "$1" == "--" ];
		then
                        shift; break
		else
                        UNRESOLVED_PARAMS="${UNRESOLVED_PARAMS} $1"
                fi
                shift
        done        
        
        FLOW_PARAMS="$@"
	
	# Try to get remaining parameters as flow parameters
	#if [ $# -ge $OPTIND ]
	#then
	#	# Shift options
	#	shift `expr $OPTIND - 1`
	#	# First parameter is .ini file
	#	INI_FILE="${1}"
	#	# Second parameter is flow parameters
	#	FLOW_PARAMS="${2}"
	#fi
	
	# Was any ini file set?
	#if [ ! -n "${INI_FILE}" ]
	#then
	#	echo ""
	#	echo "Error: No ini file"
	#	echo ""
	#	print_help
	#	exit 12
	#fi
	
	        # Was output file set?
        if [ -n "${OUT_FILE}" ]
        then
                FLOW123D_OUTPUT="${OUT_FILE}"
                # Clear output file.
                echo -n "" > "${FLOW123D_OUTPUT}"
        else
                unset FLOW123D_OUTPUT
        fi
        
}


# takes variable TIMEOUT in format HH:MM:SS
# and expands HH (hours) and MM (minutes)
function expand_timeout() {
    # echo "input TIMEOUT: $TIMEOUT"
    HH=${TIMEOUT%%:*}
    REST=${TIMEOUT#*:}
    if [ "${HH}" != "${REST}" ]
    then
      MM=${REST%%:*}
      REST=${REST#*:}
      if [ "${MM}" != "${REST}" ]
      then
        # format HH:MM:SS
        TIMEOUT=$(( 3600 * ${HH} + 60 * ${MM} + ${REST} ))
      else
        # format MM:SS
        TIMEOUT=$(( 60 * ${HH} + ${REST} ))
      fi        
    fi
    # format SS
    
    # echo "expanded TIMEOUT: $TIMEOUT"
}


###########################################################################
# Default function that run flow123d
#
# These functions can use namely variables:
#
# NP is number of procs used to compute
# MPIEXEC is relative path to bin/mpiexec
# FLOW123D is relative path to bin/flow123d (.exe)
# FLOW_PARAMS is list of parameters of flow123d
# MEM - memory limit
# QUEUE - if set, we redirect stdout and stderr to the file ${QUEUE}.<date and time>
#
# Function has to set variable ${STDOUT_FILE} with name of file 
# containing redirected stdout and stderr  of the run.
# This file has to appear after the end of computation.
#
# Function set variable EXIT_STATUS to nonzero value in the case of an error.

function call_flow() {     
      # Check if Flow123d exists and it is executable file
      if [ -x "${FLOW123D}" ]
      then
              (
              ulimit -S -v ${MEM_LIMIT}
              nice --adjustment="${NICE}" ${CALL_TIME_LIMIT_SH} "${MPIEXEC}" -np ${NP} "${FLOW123D}" ${FLOW_PARAMS}
              exit $?
              )
              return $?
      else        
              echo "Error: can not find executable: '${FLOW123D_REL}'; WORKDIR: ${WORKDIR}"
              return 11
      fi        
}


# 
function run_flow()
{
        EXIT_STATUS=0
        

	

	# Was memory limit set?
	if [ -n "${MEM}" ]
	then
		# Set up memory limits that prevent too allocate too much memory.
		# ulimit is bash commad, it accepts limit specified in kB
		MEM_LIMIT=`expr ${MEM} \* 1024`
		#echo "MEMORY LIMIT: ${MEM_LIMIT}"
	fi
	

	# form name of the redirected stdout and stderr
        if [ -n "${QUEUE}" ] 
        then 
                STDOUT_FILE="${QUEUE}.`date +%y.%m.%d_%T`" 
                rm -f `pwd`/${STDOUT_FILE}
        fi

	# Was timeout set?
	if [ -n "${TIMEOUT}" ]
	then
                expand_timeout
                CALL_TIME_LIMIT_SH="${TIME_LIMIT_SH} -t ${TIMEOUT}"
        fi    
           

        # Flow123d runs with changed priority (19 is the lowest priority)
        if [ -n "${QUEUE}" ]
        then
                call_flow  > "/tmp/${STDOUT_FILE}" 2>&1
        else
                call_flow
        fi
        EXIT_STATUS=$?
        

	if [ -n "${QUEUE}" ]
	then
                mv "/tmp/${STDOUT_FILE}" "`pwd`/${STDOUT_FILE}"
                STDOUT_FILE="`pwd`/${STDOUT_FILE}"
        fi        
	
}

# Parse command-line arguments
parse_arguments "$@"		


# If there is hostname specific script for running flow123d, then use run_flow()
# function from this script, otherwise default run_flow() function will be used
if [ -f "${0%/*}/config/${HOSTNAME}.sh" ]
then
	. "${0%/*}/config/${HOSTNAME}.sh"
fi 

# Run Flow123d
run_flow

# print the file with merged stdout, that will appear after the job is finished
if [ -n "${STDOUT_FILE}" ]; then echo "REDIRECTED: ${STDOUT_FILE}"; fi

exit ${EXIT_STATUS}
