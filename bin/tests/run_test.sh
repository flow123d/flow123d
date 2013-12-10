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
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#
# Author(s): Jiri Hnidek <jiri.hnidek@tul.cz>
#

# Syntax:
#
#       run_test.sh   "<list of input files>"  "<list of processors counts>" "parameters passed to flow" [update]
#
# 
# For every input file and every processor count run Flow123d and compare result files against saved results.
# If the parameter 'update' is given, the script do not raise an error if the result files do not match but rather
# ask user to replace reference results.
# 



#
# Note: This script depends on flow123d.sh in forcing the time out limit especially in the case that Flow123d run extremly long due to possible error.
# TODO: 
#  * allow differned setting (timeout, queue, memory limit ...) for different machines
#  * report timer from profiler (if exists)
#
#  
# set -x


# Every test has to be finished in $TIME_OUT seconds. Flow123d will be killed after
# this timeout seconds. It prevents test to run in never ending loop, when development
# version of Flow123d contains such error.
TIMEOUT=120


# Relative path to Flow123d script from the directory,
# where this script is placed
FLOW123D_SH="../flow123d.sh"
# Relative path to Flow123d binary from current/working directory
FLOW123D_SH="${0%/*}/${FLOW123D_SH}"


# Relative path to ndiff checking correctness of output files from this directory
NDIFF="../ndiff/ndiff.pl"
# Relative path to ndiff binary from current/working directory
NDIFF="${0%/*}/${NDIFF}"

# Directory that flow123d usually uses for output files
OUTPUT_DIR="./output"

# Flow output is redirected to the file. This output is printed to the output
# only in situation, when Flow123d return error (not zero) exit status
FLOW123D_OUTPUT="${OUTPUT_DIR}/flow_stdout_stderr.log"

# Ndiff output is redirected to the file. This output is printed to the output
# only in situation, when output file isn't "same" as reference output file
NDIFF_OUTPUT="./ndiff_stdout_stderr.log"

# Directory containing reference output files
REF_OUTPUT_DIR="./ref_output"

# Directory ./output containing output files of one test will be copied to the
# following directory. Each sub-directory will be named according .ini file and
# number of processes
TEST_RESULTS="./test_results"


# Variable with exit status. Possible values:
# 0 - no error, all tests were finished correctly
# 1 - some important file (flow123d, ini file) doesn't exist or permission
#     are not granted
# 2 - flow123d was not finished correctly
# 3 - execution of flow123d wasn't finished in time
EXIT_STATUS=0

# Set up memory limits (in MB) per process. Poor memory leak prevention.
# Doesn't work under Cygwin (ulimit not supported).
MEMORY_LIMIT=300


# First parameter has to be list of ini files; eg: "flow.ini flow_vtk.ini"
INI_FILES="$1"

# Second parameter has to be number of processors to run on; eg: "1 2 3 4 5"
N_PROC="$2"

# The last parameter could contain additional flow parameters
FLOW_PARAMS="$3"

# If the parameter 'update' is given, the script do not raise an error if the result files do not match but rather
# ask user to replace reference results.
UPDATE_REFERENCE_RESULTS="$4"


# set executable for awk text processor
AWK="awk"

# Following function is used to copy output files to test_results
function copy_outputs {

	# Check, if function was called with right number of arguments
	if [ $# -ge 2 ]
	then
		INI_FILE="${1}"
		NP="${2}"
	else
		echo "Error: $0 called with wrong number of arguments: $#"
		return 1
	fi

	# Remove all results of preview test
	rm -rf "${TEST_RESULTS}/${INI_FILE}.${NP}"

	# Try to create new directory for results
	mkdir -p "${TEST_RESULTS}/${INI_FILE}.${NP}" > /dev/null 2>&1

	# Was directory created successfully?
	if [ $? -eq 0 -a -d "${TEST_RESULTS}/${INI_FILE}.${NP}" ]
	then
		# Move content of ./output directory to the directory containing
		# results of tests
		cp -R "${OUTPUT_DIR}"/* "${TEST_RESULTS}/${INI_FILE}.${NP}"/
		if [ $? -ne 0 ]
		then
			echo "Error: can't copy files to: ${TEST_RESULTS}/${INI_FILE}.${NP}"
			rm -rf "${TEST_RESULTS}/${INI_FILE}.${NP}"
			return 1
		fi
	else
		echo "Error: can't create directory: ${TEST_RESULTS}/${INI_FILE}.${NP}"
		return 1
	fi
}

# Takes reference directory as $1 and current test result directory as $2.
# Updates all files and subdirs in $1 with corresponding files in $2.
# Raise an error if the source file in $2 does not exist.
function update_ref_results {
	target_dir=$1
	source_dir=$2
	#TODO: Add code that do something useful :-)
	find $target_dir -name * -exec echo '{}'
}


# Following function is used for checking one output file
function check_file {
	
	# Check, if function was called with right number of arguments
	if [ $# -ge 3 ]
	then
		INI_FILE="${1}"
		NP="${2}"
		OUT_FILE="${3}"
	else
		echo " [Failed]"
		echo "Error: $0 called with wrong number of arguments: $#"
		return 1
	fi
	
	# Print some debug information to the output of ndiff
	echo "ndiff: ${REF_OUTPUT_DIR}/${INI_FILE}/${file} ${TEST_RESULTS}/${INI_FILE}.${NP}/${OUT_FILE}" \
	>> "${TEST_RESULTS}/${INI_FILE}.${NP}/${NDIFF_OUTPUT}" 2>&1
	
	echo "" >> "${TEST_RESULTS}/${INI_FILE}.${NP}/${NDIFF_OUTPUT}" 2>&1
	
	# Compare output file using ndiff
	${NDIFF} \
		"${REF_OUTPUT_DIR}/${INI_FILE}/${OUT_FILE}" \
		"${TEST_RESULTS}/${INI_FILE}.${NP}/${OUT_FILE}" \
		>> "${TEST_RESULTS}/${INI_FILE}.${NP}/${NDIFF_OUTPUT}" 2>&1
	
	# Check result of ndiff
	if [ $? -eq 0 ]
	then
		echo -n "."
		return 0
	else
		echo " [Failed]"
		echo "Error: file ${TEST_RESULTS}/${INI_FILE}.${NP}/${OUT_FILE} is too different."
		return 1
	fi
}

# Following function is used for checking output files
function check_outputs {

	echo -n "Checking output files ."

	# Check, if function was called with right number of arguments
	if [ $# -ge 2 ]
	then
		INI_FILE="${1}"
		NP="${2}"
	else
		echo " [Failed]"
		echo "Error: $0 called with wrong number of arguments: $#"
		return 1
	fi

	
	# Does exist reference directory for this .ini file?
	if [ ! -d "${REF_OUTPUT_DIR}/${INI_FILE}" ]
	then
		echo " [Failed]"
		echo "Error: directory with reference output files doesn't exist: ${REF_OUTPUT_DIR}/${INI_FILE}"
		return 1
	fi

	REF_FILE_NUM=`ls "${REF_OUTPUT_DIR}/${INI_FILE}" | wc -l`
	# Does reference directory contain any reference file?
	if [ ${REF_FILE_NUM} -eq 0 ]
	then
		echo " [Failed]"
		echo "Error: directory ${REF_OUTPUT_DIR}/${INI_FILE} doesn't contain any reference output files."
		return 1
	fi
	unset REF_FILE_NUM

	# Does exist output directory for this .ini file and number of processes?
	if ! [ -d "${TEST_RESULTS}/${INI_FILE}.${NP}" ]
	then
		echo " [Failed]"
		echo "Error: directory with output files doesn't exist: ${TEST_RESULTS}/${INI_FILE}.${NP}"
		return 1
	fi
	
	# If an update should be performed 
	NEED_UPDATE=

	# For every file in reference directory try to find generated file
	# and do ndiff
	for file in `ls "${REF_OUTPUT_DIR}/${INI_FILE}/"`
	do
		# Required output file exists?
		if [ -e "${TEST_RESULTS}/${INI_FILE}.${NP}/${file}" ]
		then
			# Is test file directory?
			if [ -d "${TEST_RESULTS}/${INI_FILE}.${NP}/${file}" ]
			then
				# If $file is directory, then check all files in this directory
				subdir="${file}"
				for subfile in `ls "${REF_OUTPUT_DIR}/${INI_FILE}/${subdir}/"`
				do
					if [ -e "${TEST_RESULTS}/${INI_FILE}.${NP}/${subdir}/${subfile}" ]
					then
						# Compare file with reference
						check_file "${INI_FILE}" "${NP}" "${subdir}/${subfile}"
						if [ $? -ne 0 ]
						then
							if [ "${UPDATE_REFERENCE_RESULTS}" = "update" ]
							then
								NEED_UPDATE=1
							else
								return 1
							fi
						fi
					else
						echo " [Failed]"
						echo "Error: file ${TEST_RESULTS}/${INI_FILE}.${NP}/${subdir}/${subfile} doesn't exist"
						return 1
					fi
				done
			else
				# Compare file with reference
				check_file ${INI_FILE} ${NP} ${file}
				if [ $? -ne 0 ]
				then
					if [ "${UPDATE_REFERENCE_RESULTS}" = "update" ]
					then
						NEED_UPDATE=1
					else
						return 1
					fi
				fi
			fi
		else
			echo " [Failed]"
			echo "Error: file ${TEST_RESULTS}/${INI_FILE}.${NP}/${file} doesn't exist"
			return 1
		fi
	done
	
	if [ -n "${NEED_UPDATE}" ]
	then
		echo "Do you want to update reference results? (y/n) [no]"
		read UPDATE
		if [ "UPDATE" = "y" ] 
		then
			update_ref_results ${REF_OUTPUT_DIR}/${INI_FILE} ${TEST_RESULTS}/${INI_FILE}.${NP}
		fi
	fi
}


function wait_for_flow_script {
        # Wait for (finished) flow script and get its exit status ('wait' returns status of the sub process)   
        if [ -z "${FLOW_EXIT_STATUS}" ]
        then
          wait ${FLOW123D_PID}
          FLOW_EXIT_STATUS=$?
        fi  

        # set up, that flow123d was finished in time
        STDOUT_FILE=`cat "${FLOW_SCRIPT_STDOUT}" | grep "REDIRECTED: "`
        STDOUT_FILE="${STDOUT_FILE#REDIRECTED: }"  
        #echo "Waiting for ${STDOUT_FILE}."
}

#########################################################################################################################333
# MAIN


# Check if Flow123d exists and it is executable file
if ! [ -x "${FLOW123D_SH}" ]
then
	echo "Error: can't execute ${FLOW123D_SH}"
	exit 1
fi

# Make output directory
if [ ! -d ${OUTPUT_DIR} ]
then 
  mkdir -p ${OUTPUT_DIR}
fi  


FLOW_SCRIPT_STDOUT="`pwd`/flow_script.stdout"

# For every ini file run one test
for INI_FILE in $INI_FILES
do
	# Check if it is possible to read ini file
	if ! [ -e "${INI_FILE}" -a -r "${INI_FILE}" ]
	then
		echo "Error: can't read ${INI_FILE}"
		EXIT_STATUS=1
		continue 1
	fi

	for NP in ${N_PROC}
	do

		# Erase content of ./output directory
		rm -rf "${OUTPUT_DIR}"/*

		# Reset timer
		TIMER="0"

		# Flow123d runs with changed priority (19 is the lowest priority)
                "${FLOW123D_SH}" --nice 10 --mem ${MEMORY_LIMIT} -np ${NP} -ppn 1 --walltime ${TIMEOUT} -q "short" -- -s "${INI_FILE}" ${FLOW_PARAMS} >"${FLOW_SCRIPT_STDOUT}" &
                # Get PID 
		FLOW123D_PID=$!

		echo -n "Running flow123d [proc:${NP}] ${INI_FILE} ."
		IS_RUNNING=1

		# Wait max TIMEOUT seconds, then flow123d processes should be killed
		while [ ${TIMER} -lt ${TIMEOUT} ]
		do
			TIMER=`expr ${TIMER} + 1`
			echo -n "."
			#ps -o "%P %p"
			sleep 1

			# Is flow script still running?
			if [  ${IS_RUNNING} -eq 1 ]
			then 
                              ps | ${AWK} '{ print $1 }' | grep -q "${FLOW123D_PID}"
                              if [ $? -ne 0 ]
                              then
                                      IS_RUNNING="2"
                                      wait_for_flow_script
                              fi                              
			else
                              # wait for flow to finish 
                              if [ -e "${STDOUT_FILE}" ]
                              then
                                      IS_RUNNING="0"
                                      break
                              fi
			fi
		done
		
		# we wait in order to get STDOUT_FILE even in case of error (at least for interactive runs)
		wait_for_flow_script
		if [ -e ${STDOUT_FILE} ]
		then
                    mv "${STDOUT_FILE}" "${FLOW123D_OUTPUT}"
                else
                    FLOW_EXIT_STATUS=100
                fi
		
		rm -f ${FLOW_SCRIPT_STDOUT}

		# In all cases copy content of ./output to ./test_results directory
		copy_outputs "${INI_FILE}" "${NP}"

		

		# Was Flow123d finished correctly?
		if [ ${IS_RUNNING} -eq 0 -a ${FLOW_EXIT_STATUS} -eq 0 ]
		then
			echo " [Success:${TIMER}s]"
			
			# Check correctness of output files
			check_outputs "${INI_FILE}" "${NP}"

			# Were all output files correct?
			if [ $? -eq 0 ]
			then
				echo " [Success]"
			else
				EXIT_STATUS=10
				# Try next ini file
				continue 2 
			fi
		else
			echo " [Failed:error]"
			EXIT_STATUS=1
			# No other test will be executed
			break 2
		fi
	done
done

# Print redirected stdout to stdout only in situation, when some error occurred
if [ $EXIT_STATUS -gt 0 -a $EXIT_STATUS -lt 10 ]
then
	echo "Error in execution: ${FLOW123D_SH} -s ${INI_FILE} ${FLOW_PARAMS}"
	cat "${FLOW123D_OUTPUT}"
fi

rm -f "${FLOW123D_OUTPUT}"

exit ${EXIT_STATUS}

