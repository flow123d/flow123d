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
# Author(s): Jiri Hnidek <jiri.hnidek@tul.cz>
#


#
# Note: This script assumes that Flow123d can contain any error. It means that
# there could be nevere ending loop, flow could try to allocate infinity
# amount of memory, etc.
#


# Every test has to be finished in 60 seconds. Flow123d will be killed after
# 60 seconds. It prevents test to run in never ending loop, when developemnt
# version of Flow123d contains such error.
TIMEOUT=120

# Try to use MPI environment variable for timeout too. Some implementation
# of MPI supports it and some implementations doesn't.
export MPIEXEC_TIMEOUT=${TIMEOUT}

# Flow output is redirected to the file. This output is printed to the output
# only in situation, when Flow123d return error (not zero) exit status
FLOW123D_OUTPUT="./flow_stdout.log"

# Relative path to Flow123d binary from the directory,
# where this script is placed
FLOW123D="../flow123d"
# Relative path to Flow123d binary from current/working directory
FLOW123D="${0%/*}/${FLOW123D}"

# Relative path to mpiexec from the directory, where this script is placed
MPIEXEC="../mpiexec"
# Relative path to mpiexec binary from current/working directory
MPIEXEC="${0%/*}/${MPIEXEC}"

# Relative path to ndiff checking corectness of output files from this directory
NDIFF="../ndiff/ndiff.pl"
# Relative path to ndiff binary from current/working directory
NDIFF="${0%/*}/${NDIFF}"


# Variable with exit status. Possible values:
# 0 - no error, all tests were finished correctly
# 1 - some important file (flow123d, ini file) doesn't exist or presmission
#     are not granted
# 2 - flow123d was not finished corectly
# 3 - execution of flow123d wasn't finished in time
EXIT_STATUS=0

# Set up memory limits that prevent too allocate too much memory.
# The limit of virtual memory is 200MB (memory, that could be allocated)
ulimit -S -v 200000


# First parameter has to be list of ini files; eg: "flow.ini flow_vtk.ini"
INI_FILES="$1"

# Secons parameter has to be number of processors to run on; eg: "1 2 3 4 5"
N_PROC="$2"

# The last parameter could contain additional flow params
FLOW_PARAMS="$3"


# Following function is used for checking output files
function check_output {
	${NDIFF} ${1} ${2}
}


# Check if Flow123d exists and it is executable file
if ! [ -x "${FLOW123D}" ]
then
	echo "Error: can't execute ${FLOW123D}"
	exit 1
fi

# Check if mpiexec exists and it is executable file
if ! [ -x "${MPIEXEC}" ]
then
	echo "Error: can't execute ${MPIEXEC}"
	exit 1
fi

# For every ini file run one test
for INI_FILE in $INI_FILES
do
	# Check if it is possible to read ini file
	if ! [ -e "${INI_FILE}" -a -r "${INI_FILE}" ]
	then
		echo "Error: can't read ${INI_FILE}"
		exit 1
	fi

	for NP in ${N_PROC}
	do
		# Clear output file for every new test. Output of passed test isn't
		# important. It is usefull to see the output of last test that failed.
		echo "" > "${FLOW123D_OUTPUT}"

		# Reset timer
		TIMER="0"

		# Flow123d runs with changed priority (19 is the lowest priority)
		nice --adjustment=10 "${MPIEXEC}" -np ${NP} "${FLOW123D}" -S "${INI_FILE}" ${FLOW_PARAMS} > "${FLOW123D_OUTPUT}" 2>&1 &
		PARENT_MPIEXEC_PID=$!
		# Fait for child proccesses
		sleep 1
		MPIEXEC_PID=`ps -o "${PARENT_MPIEXEC_PID} %P %p" | gawk '{if($1==$2){ print $3 };}'`

		echo -n "Runing flow123d [proc:${NP}] ${INI_FILE} ."
		IS_RUNNING=1

		# Wait max TIMEOUT seconds, then kill mpiexec, that executed flow123d proccesses
		while [ ${TIMER} -lt ${TIMEOUT} ]
		do
			TIMER=`expr ${TIMER} + 1`
			echo -n "."
			#ps -o "%P %p"
			sleep 1

			# Is mpiexec and flow123d still running?
			ps | gawk '{ print $1 }' | grep -q "${PARENT_MPIEXEC_PID}"
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
			echo " [Failed:loop]"
			# Send SIGTERM to mpiexec. It should kill child proccessed
			kill -s SIGTERM ${MPIEXEC_PID} #> /dev/null 2>&1
			EXIT_STATUS=2
			# No other test will be executed
			break 2
		else
			# Get exit status variable of mpiexec executing mpiexec executing flow123d
			wait ${PARENT_MPIEXEC_PID}
			MPIEXEC_EXIT_STATUS=$?

			# Was Flow123d finished corectly?
			if [ ${MPIEXEC_EXIT_STATUS} -eq 0 ]
			then
				echo " [Success:${TIMER}s]"
			else
				echo " [Failed:error]"
				EXIT_STATUS=1
				# No other test will be executed
				break 2
			fi
		fi
	done
done

# Print redirected stdout to stdout only in situation, when some error ocured
if [ $EXIT_STATUS -ne 0 ]
then
	echo "Error in execution: ${FLOW123D} -S ${INI_FILE} -- ${FLOW_PARAMS}"
	cat "${FLOW123D_OUTPUT}"
fi

rm -f "${FLOW123D_OUTPUT}"

exit ${EXIT_STATUS}

