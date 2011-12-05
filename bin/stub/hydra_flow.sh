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
# $Id: hydra_make_pbs.sh 834 2011-01-02 15:54:35Z michal.nekvasil $
# $Revision: 834 $
# $LastChangedBy: michal.nekvasil $
# $LastChangedDate: 2011-01-02 16:54:35 +0100 (ne, 02 I 2011) $
#

# Warning: This script is stub. Do not use it!

# NPROC is defined in run_flow.sh and its number of procs used to compute
# MPI_RUN is defined in run_flow.sh and its relative path to bin/mpiexec
# EXECUTABLE is defined in run_flow.sh and its relative path to bin/flow123d (.exe)
# FLOW_PARAMS is defined in run_flow.sh
# INI is name of inifile and its defined in run_flow.sh

# Some important files
export LOCK_FILE="./lock"
export ERR_FILE="err.log"
export OUT_FILE="out.log"

# Copy following text to the file hydra_flow.qsub
# ===============================================
cat << xxEOFxx > hydra_flow.qsub
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

# What is it!?
export OMPI_MCA_plm_rsh_disable_qrsh=1

# Execute Flow123d using mpirun
"$MPI_RUN" -np $NPROC "$EXECUTABLE" -S "$INI" $FLOW_PARAMS 2>${ERR_FILE} 1>${OUT_FILE}

# Delete lock file
rm -f ${LOCK_FILE}

# End of hydra_flow.qsub
xxEOFxx
# ===============================================

# Test, if lock file exist
if [ -e ${LOCK_FILE} ] 
then
	echo "Error: the working directory is locked by onother PBS job!"
	exit 1
else
	# Create new lock file
	touch ${LOCK_FILE}
fi

# Add new PBS job to the queue
qsub -pe orte $NPROC hydra_flow.qsub
