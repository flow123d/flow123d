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

set -x 

# Relative path to mpiexec from the directory, where this script is placed
MPIEXEC="./mpiexec"
# Relative path to mpiexec binary from current/working directory
MPIEXEC="${0%/*}/${MPIEXEC}"

# Relative path to Flow123d binary from the directory,
# where this script is placed
FLOW123D="./flow123d"
# Absolute path to Flow123d binary from current/working directory
FLOW123D="${PWD}/${0%/*}/${FLOW123D}"
# Absolute path to working directory
WORKDIR="${PWD}"

AWK="awk"

# Print help to this script
function print_help {
	echo "SYNTAX: flow123d.sh [OPTIONS] INI_FILE [\"FLOW_PARAMS\"]"
	echo ""
	echo "INI_FILE            Flow123d will load configuration from INI_FILE"
	echo "FLOW_PARAMS         Flow123d will use it's specific parameters."
	echo "                    Such parameters has to be in brackets"
	echo ""
	echo "OPTIONS:"
	echo "    -h              Print this help"
	echo "    -t TIMEOUT      Flow123d can be executed only TIMEOUT seconds"
	echo "    -m MEM          Flow123d can use only MEM bytes"
	echo "    -n NICE         Run Flow123d with changed (lower) priority"
	echo "    -p N            Run Flow123d using N parallel procces" 
	echo "    -r OUT_FILE     Stdout and Stderr will be redirected to OUT_FILE"
	echo "    -s              Working directory will be current directory (default)"
	echo "    -S              Working directory will be relative path to ini file"
	echo "    -q QUEUE        Name of queue to use for batch processing (if supported by system)"
	echo ""
	echo "RETURN VALUES:"
	echo "    0               Flow123d process exited normaly"
	echo "    10              Bad argument/option"
	echo "    11              Flow123d can't be executed"
	echo "    12              INI file doesn't exist"
	echo "    13              Flow123d can't read INI file"
	echo "    14              Flow123d process was killed (TIMEOUT exceeded)"
	echo "    15              Flow123d process crashed"
	echo ""
}

# Parse script arguments
function parse_arguments()
{
	# Make sure that INI_FILE is not set
	unset INI_FILE

	# Default number of processes
	NP=1
	
	# Default behavior is to use -s parameter
	FLOW_OPT="-s"
	
	# Parse arguments with bash builtin command getopts
	while getopts ":ht:m:n:p:q:r:sS" opt
	do
		case ${opt} in
		h)
			print_help
			exit 0
			;;
		t)
			TIMEOUT="${OPTARG}"
			;;
		m)
			MEM="${OPTARG}"
			;;
		n)
			NICE="${OPTARG}"
			;;
		p)
			NP="${OPTARG}"
			;;
                q)
                        QUEUE="${OPTARG}"
                        ;;
		r)
			OUT_FILE="${OPTARG}"
			;;
		s)
			FLOW_OPT="-s"
			;;
		S)
			FLOW_OPT="-S"
			;;
		\?)
			echo ""
			echo "Error: Invalid option: -$OPTARG"
			echo ""
			print_help
			exit 10
			;;
		esac
	done
	
	# Try to get remaining parameters as flow parameters
	if [ $# -ge $OPTIND ]
	then
		# Shift options
		shift `expr $OPTIND - 1`
		# First parameter is .ini file
		INI_FILE="${1}"
		# Second parameter is flow parameters
		FLOW_PARAMS="${2}"
	fi
	
	# Was any ini file set?
	if [ ! -n "${INI_FILE}" ]
	then
		echo ""
		echo "Error: No ini file"
		echo ""
		print_help
		exit 12
	fi
}

# Default function that run flow123d
function run_flow()
{
	# Check if Flow123d exists and it is executable file
	if ! [ -x "${FLOW123D}" ]
	then
		echo "Error: can't execute ${FLOW123D}"
		exit 11
	fi
	
	# Check if it is possible to read ini file
	if ! [ -e "${INI_FILE}" -a -r "${INI_FILE}" ]
	then
		echo "Error: can't read ${INI_FILE}"
		exit 13
	fi
	
	# Was memory limit set?
	if [ -n "${MEM}" ]
	then
		# Set up memory limits that prevent too allocate too much memory.
		# The limit of virtual memory is 200MB (memory, that could be allocated)
		MEM_LIMIT=`expr ${MEM} \* 1000`
		ulimit -S -v ${MEM_LIMIT}
	fi
	
	# Was nice set?
	if [ ! -n "${NICE}" ]
	then
		NICE=0
	fi
	
	# Was output file set?
	if [ -n "${OUT_FILE}" ]
	then
		FLOW123D_OUTPUT="${OUT_FILE}"
		# Clear output file.
		echo -n "" > "${FLOW123D_OUTPUT}"
	else
		unset FLOW123D_OUTPUT
	fi
	
	# Was timeout set?
	if [ -n "${TIMEOUT}" ]
	then
		# Flow123d runs with changed priority (19 is the lowest priority)
		if [ -n "${FLOW123D_OUTPUT}" ]
		then
			nice --adjustment="${NICE}" "${FLOW123D}" ${FLOW_OPT} "${INI_FILE}" ${FLOW_PARAMS} > "${FLOW123D_OUTPUT}" 2>&1 &
		else
			nice --adjustment="${NICE}" "${FLOW123D}" ${FLOW_OPT} "${INI_FILE}" ${FLOW_PARAMS} &
		fi
		FLOW123D_PID=$!
		
		IS_RUNNING=1

		TIMER=0
		# Wait max TIMEOUT seconds, then kill flow123d processes
		while [ ${TIMER} -lt ${TIMEOUT} ]
		do
			TIMER=`expr ${TIMER} + 1`
			sleep 1

			# Is flow123d process still running?
			ps | ${AWK} '{ print $1 }' | grep -q "${FLOW123D_PID}"
			if [ $? -ne 0 ]
			then
				# set up, that flow123d was finished in time
				IS_RUNNING="0"
				break 1
			fi
		done
		
		# Was Flow123d finished during TIMEOUT or is it still running?
		if [ ${IS_RUNNING} -eq 1 ]
		then
			# Send SIGTERM to flow123d.
			kill -s SIGTERM ${FLOW123D_PID} #> /dev/null 2>&1
			exit 14
		else
			# Get exit status variable of flow123d
			wait ${FLOW123D_PID}
			FLOW123D_EXIT_STATUS=$?

			# Was Flow123d finished correctly?
			if [ ${FLOW123D_EXIT_STATUS} -eq 0 ]
			then
				exit 0
			else
				exit 15
			fi
		fi
	else
		if [ -n "${FLOW123D_OUTPUT}" ]
		then
			nice --adjustment="${NICE}" "${FLOW123D}" ${FLOW_OPT} "${INI_FILE}" ${FLOW_PARAMS} > "${FLOW123D_OUTPUT}" 2>&1
		else
			nice --adjustment="${NICE}" "${FLOW123D}" ${FLOW_OPT} "${INI_FILE}" ${FLOW_PARAMS}
		fi
	fi
	
	exit 0
}

# Parse command-line arguments
parse_arguments "$@"		

# If there is hostname specific script for running flow123d, then use run_flow()
# function from this script, otherwise default run_flow() function will be used
if [ -f "${0%/*}/config/${HOSTNAME//./_}.sh" ]
then
	. "${0%/*}/config/${HOSTNAME//./_}.sh"
fi 

# Run Flow123d
run_flow
