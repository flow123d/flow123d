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

# Relative path to Flow123d binary from the directory,
# where this script is placed
FLOW123D="./flow123d"
# Relative path to Flow123d binary from current/working directory
FLOW123D="${0%/*}/${FLOW123D}"


# Print help to this script
function print_help {
	echo "SYNTAX: flow123d.sh [OPTIONS] -i INI_FILE [ -f FLOW_PARAMS]"
	echo ""
	echo "OPTIONS:"
	echo "    -h              Print this help"
	echo "    -t TIMEOUT      Flow123d can be executed only TIMEOU seconds"
	echo "    -m MEM          Flow123d can use only MEM bytes"
	echo "    -i INI_FILE     Flow123d will load configuration from INI_FILE"
	echo "    -f FLOW_PARAMS  Flow123d will use it's specific params"
	echo ""
}

# Make sure that INI_FILE is not set
unset INI_FILE

# Parse arguments with bash builtin command getopts
while getopts ":ht:m:i:f:" opt
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
	i)
		# -i has to be followed by ini file; e.g.: "flow.ini"
		INI_FILE="${OPTARG}"
		;;
	f)
		# Additional flow parametres; e.g.: "-i input -o output"
		FLOW_PARAMS="${OPTARG}"
		;;
	\?)
		echo ""
		echo "Error: Invalid option: -$OPTARG"
		echo ""
		print_help
		exit 1
		;;
	esac
done

# Was any ini file set?
if [ ! -n "${INI_FILE}" ]
then
	echo ""
	echo "Error: no ini file"
	echo ""
	print_help
	exit 1
else
	#"${FLOW123D}" -S "${INI_FILE}" ${FLOW_PARAMS}
	echo "timeout: ${TIMEOUT} mem_limit: ${MEM} ini_file: ${INI_FILE} ${FLOW_PARAMS}"
fi


