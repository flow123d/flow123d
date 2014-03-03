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

# Function that is used for running flow123d at hydra cluster
function run_flow()
{
	# Some important files
	export ERR_FILE="err.log"
	export OUT_FILE="out.log"
	
	QSUB_FILE="/tmp/${USER}-flow123.qsub"
	if [ -z "${QUEUE}" ]; then QUEUE=normal; fi

	rm -f ${QSUB_FILE}
			
# Copy following text to the file /tmp/firstname.surname-hydra_flow.qsub
# ======================================================================
cat << xxEOFxx > ${QSUB_FILE}
#!/bin/bash
#
# Specific PBS setting
#
#PBS -S /bin/bash 
#PBS -j oe 
#PBS -N flow123d
#PBS -m bae
##################BS -l walltime=24:00:00
#PBS -l select=1:ncpus=$NP:host=rex
#PBS -l place=free:shared 

# set paths to intel libraries
source /opt/intel/Compiler/11.1/046/bin/iccvars.sh ia64 

# set same LD paths as in the evironment from which the task is started
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ARCH_LD_LIBRARY_PATH

# load module for SGI MPI library
. /usr/share/modules/init/bash
module load mpt

# ???
export KMP_MONITOR_STACKSIZE=64K   

cd ${WORKDIR}

# header - info about task
uname -a
echo JOB START: `date` 
pwd

# dplace is used for better process placement
dplace mpirun -np $NP "$FLOW123D" $FLOW_PARAMS 
# 2>${ERR_FILE} 1>${OUT_FILE}
	
# End of flow123d.qsub
xxEOFxx
# ======================================================================

	if [ -f ${QSUB_FILE} ]
	then    
		# Add new PBS job to the queue
		echo "qsub -r n -p 1023 -q ${QUEUE} ${QSUB_FILE}"
		qsub -r n -p 1023 -q ${QUEUE} ${QSUB_FILE}
		# Remove obsolete script
		rm ${QSUB_FILE}
	else
		exit 1
	fi		
	
	exit 0
}