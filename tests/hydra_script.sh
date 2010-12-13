#!/bin/bash

echo "
#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
 
export OMPI_MCA_plm_rsh_disable_qrsh=1
$HOME/flow/1.6.0_modular/bin/mpiexec $NSLOTS $HOME/flow/1.6.0_modular/bin/$EXECUTABLE -s $INI 2>err 1>out" >run_hydra.qsub

chmod u+x hydra.qsub