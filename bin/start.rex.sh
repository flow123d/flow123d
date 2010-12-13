#!/bin/bash
#
# Starting script for sgi REX, CIV, Prague
# 
# For help run:
# start.rex.sh --help

set -x 

# direcotry of starting script - s3d binary is relative to this path
SCRIPT_PATH="`pwd`/${0%/*}"
#CURRENT_SCRIPT="`pwd`/$0
source $SCRIPT_PATH/start.common

ARCHHELP="
  -q [queue]  start job noninteractively in 'queue' of the batch system
	      (on Rex omit queue name)
  
  -wt time    set upper esstimate of the time your job need to finish,
              is necessary for the automatic choice of a queue
              default value is 4:00:00 (i.e. 4 hours)
                
  -m          set MPI_NAP to allow run multiple processes on one processor
              (slow workaround some bug in SGI's MPI)

  \$PBS_USR_PARAMS  - by this variable you can specify set of parameters
              you wish to be passed to qsub every time
              (usefull for -m and -M parametr for set mail reporting
               see: man qsub)
"

NOHUP=""
# wall time - PBS use this to choose the queue
WTIME="4:00:00"
QUEUE=
# here you can spacify additional dirs of shared libraries
#ARCH_LD_LIBRARY_PATH=/opt/intel/mkl/8.0.1/lib/64
ARCH_LD_LIBRARY_PATH=

SCREEN_REDIRECT_FILE=run.out 	# should contain date


# parse 'rex' specific parameters
ParseArch() {
while [ \( -n "$1" \) -a \( ! "$1" == "--" \) ]
do
  if [ "$1" == "-wt" ]; then
    shift
    if [ \( -n "$1" \) -a \( "${1#-}" == "${1}" \) ]
    then 
      WTIME=$1
      shift
    else
      echo "Unspecified wall time ... ignore"
    fi
  elif [ "$1" == "-q" ]; then
    shift
    NOHUP="yes"
#    if [ \( -n "$1" \) -a \( "${1#-}" == "${1}" \) ]
#    then 
#      QUEUE=$1
#      shift
#    else
#      echo "Unspecified queue."
#      exit 1
#    fi
  elif [ "$1" == "-m" ]; then
    export MPI_NAP=1
    shift
  else
    shift
  fi
done  
}


# make script s3dpar.job 
# this can be used to start s3d over PBS or directly
# this script has no parameters !!
MakeStartScript () {
echo "
#!/bin/bash
#
# Specific PBS setting
#
#PBS -S /bin/bash 
#PBS -j oe 
#PBS -N s3dpar
#PBS -m bae
#PBS -l walltime=$WTIME
#PBS -l select=1:ncpus=$NPROC:host=rex
#PBS -l place=free:shared 
. /usr/share/modules/init/bash
module load mpt
export KMP_MONITOR_STACKSIZE=64K   

cd $WRK_DIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ARCH_LD_LIBRARY_PATH
uname -a
echo JOB START: `date` 
pwd
mpirun -np $NPROC $SCRIPT_PATH/s3dpar $SOLVEROPT -s3d_input $START_INDIR -s3d_output $START_OUTDIR > $START_OUTDIR/$SCREEN_REDIRECT_FILE &&
$SCRIPT_PATH/outmerge $START_OUTDIR
" >s3dpar.job
chmod u+x s3dpar.job
}

RunJob() {
  set -x
  
  MakeStartScript
  if [ "$NOHUP" == "yes" ]
  then 
    # not rerunnable job; highest priority
    qsub -r n -p 1023 $PBS_USR_PARAMS ./s3dpar.job
  else
    ./s3dpar.job
  fi
}

# start parallel job
StartArchPar() {
    RunJob
}  

# start sequentional job
StartArchSeq() {
    RunJob
}

# ================================================================= MAIN

ParseCommonArgs $@
echo $NPROC
ParseArch $ARGREST
PrepareIO		# call restartor save old output

if [ $NPROC -gt 1 ] 
then
  StartArchPar
else
  StartArchSeq
fi
