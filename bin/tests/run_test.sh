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
# This script assumes that Flow123d can contain any error. It means, that
# there could be nevere ending loop, flow could try to allocate infinity
# amount of memory, etc.
#
# This script has to be run from root of source files.




# Every test has to be finished in 60 seconds. Flow123d will be killed after
# 60 seconds. It prevents test to run in never ending loop, when developemnt
# version of Flow123d contains such error.
TIMEOUT=60

# Variable with exit status. Possible values: 0 - no error; 1 - flow123d was
# not finished corectly; 2 - execution of flow123d wasn't finished in time
EXIT_STATUS=0

# First parameter has to be list of ini files; eg: "flow.ini flow_vtk.ini"
INI_FILES="$1"

# Secons parameter has to be number of processors to run on; eg: "1 2 3 4 5"
NPROC="$2"

# The last parameter could contain additional flow params
FLOW_PARAMS="$3"

# Flow123d binary
FLOW123D=./bin/flow123d




# For every ini file run one test
for INI_FILE in $INI_FILES
do
	echo "Runing flow123 ..."
	$FLOW123D -S "$INI_FILE" -- "$FLOW_PARAMS" > output_pokus.txt 2>&1 &
	FLOW123D_PID=$!
	IS_RUNNING=1

	# Wait max TIMEOUT seconds, then kill Flow123d
	while [ ${TIMEOUT} -gt 0 ]
	do
		TIMEOUT=`expr ${TIMEOUT} - 1`
		echo -n "."
		sleep 1

		# Is Flow123d still running?
		ps | gawk '{ print $1 }' | grep -q "${FLOW123D_PID}"
		if [ $? -ne 0 ]
		then
			# Flow123d was finished in time
			IS_RUNNING="0"
			break
		fi
	done

	# Was RUNNER finished during TIMEOUT or is it still running?
	if [ ${RUNNING} -eq 1 ]
	then
		echo " [Failed]"
		kill -9 ${RUNNER_PID}
		EXIT_STATUS=2
		# No other test will be executed
		break
	else
		# Get exit status variable of Flow123dd
		wait ${RUNNER_PID}

		# Was Flow123d finished corectly?
		if [ $? -eq 0 ]
		then
			echo " [Success]"
		else
			echo " [Failed]"
			EXIT_STATUS=1
			# No other test will be executed
			break
		fi
	fi
done

exit ${EXIT_STATUS}

