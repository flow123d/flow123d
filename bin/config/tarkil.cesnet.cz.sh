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
# $Id: hydra_kai_tul_cz.sh 1527 2012-02-08 07:30:20Z jiri.hnidek $
# $Revision: 1527 $
# $LastChangedBy: jiri.hnidek $
# $LastChangedDate: 2011-01-02 16:54:35 +0100 (ne, 02 I 2011) $
#

# NP is number of procs used to compute
# MPIEXEC is relative path to bin/mpiexec
# FLOW123D is relative path to bin/flow123d (.exe)
# FLOW_PARAMS is list of parameters of flow123d
# INI_FILE is name of .ini file
# WORKDIR is directory from which flow123d.sh was started
# TIMEOUT is max time to run
#
# sets variable STDOUT_FILE to the file name of the joined redirected stdout and stderr
#
# TODO: 
# * the job is started from /storage/.../home/USER/...
#   It could be worth of trying to copy i/o data to/from scratch disk; it may be faster for seek operations.
# * this script depends on build with particular combination of modules and on the debian 6.0 system 
#   one possible enhancement: try to compile on debian 5.0 
#
# set -x 

# Function that is used for running flow123d at hydra cluster
function run_flow()
{
        EXIT_STATUS=0
        
	# Some important files
	export ERR_FILE="err.log"
	export OUT_FILE="out.log"
	
	QSUB_FILE="/tmp/${USER}-flow123.qsub"
        rm -f ${QSUB_FILE}
	
	if [ -z "${QUEUE}" ]; then QUEUE=normal; fi
	if [ -z "${PPN}" ]; then PPN=2; fi
        if [ -z "${MEM}" ]; then MEM="$(( ${PPN} * 2))"; fi
        # divide and round up
        NNodes="$(( ( ${NP} + ${PPN} -1 ) / ${PPN} ))"
        if [ -n "${TIMEOUT}" ]; then SET_WALLTIME="-l walltime=${TIMEOUT}";fi
			
# Copy following text to the file /tmp/firstname.surname-hydra_flow.qsub
# ======================================================================
cat << xxEOFxx > ${QSUB_FILE}
#!/bin/bash
#
# Specific PBS setting
#
#PBS -S /bin/bash 
#PBS -N flow123d
#PBS -j oe
#################
# load modules
module add svn-1.7.6
module add intelcdk-12
module add boost-1.49
module add mpich-p4-intel
module add cmake-2.8
module add python-2.6.2
module unload mpiexec-0.84
module unload mpich-p4-intel
module add openmpi-1.6-intel



cd ${WORKDIR}

# header - info about task
uname -a
echo JOB START: `date` 
pwd

echo mpirun -n $NP "$FLOW123D" $FLOW_PARAMS  

mpirun -n $NP "$FLOW123D" $FLOW_PARAMS  
  
	
# End of flow123d.qsub
xxEOFxx
# ======================================================================

	if [ -f ${QSUB_FILE} ]
	then    
                OPTIONS="-l nodes=${NNodes}:ppn=${PPN}:x86_64:nfs4:debian60 -l mem=${MEM}mb ${SET_WALLTIME} ${UNRESOLVED_PARAMS} -q ${QUEUE}"
		# Add new PBS job to the queue
		echo "qsub ${OPTIONS} ${QSUB_FILE}"
		
		JOB_NAME=`qsub ${OPTIONS} ${QSUB_FILE}`
		if [ $? -ne 0 ]; then return 1;fi
		
                # construct STDOUT_FILE
                JOB_NAME=${JOB_NAME%%.*}
                echo "job number: ${JOB_NAME}"
                STDOUT_FILE="flow123d.o${JOB_NAME}"
		# Remove obsolete script
		# rm ${QSUB_FILE}
	else
                # can not create QSUB script
                EXIT_STATUS=1
	fi		
	
}
