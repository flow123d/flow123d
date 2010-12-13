#!/bin/bash

#SCRIPT_DIR_PATH is defined in run_flow.sh and it's absolut path to dir wehere is flow_run.sh
#MPI_RUN is defined in run_flow.sh and its relative path to bin/mpiexec
#EXECUTABLE is defined in run_flow.sh and its relative path to bin/flow123d (.exe)

echo "
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
 
cd $SCRIPT_DIR_PATH
touch lock
export OMPI_MCA_plm_rsh_disable_qrsh=1
$SCRIPT_DIR_PATH/$MPI_RUN $NSLOTS $SCRIPT_DIR_PATH/$EXECUTABLE -s $INI 2>err 1>out
rm lock" >hydra_run_pbs.qsub

chmod u+x hydra_run_pbs.qsub