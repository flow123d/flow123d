#!/bin/bash

if [ "$1" == "-t" ]
then
  shift
  TIME=$1
  shift
else
  echo "Usage: time_limit.sh -t <time> <command>"
  echo "Run <command> and kill it after <time> seconds."
fi

# run command on background, take its PID
#echo "$@"
$@ &
COMMAND_PID=$!

# wait until COMMAND is finished or timeout
END_TICK="${TIME}0"
while [ -e /proc/${COMMAND_PID} -a "${TICK}" != "${END_TICK}" ]
do 
  sleep 0.1
  TICK=$(($TICK+1))
done

# kill possibly running COMMAND
if [ -e /proc/${COMMAND_PID} ]
then
  kill ${COMMAND_PID}
  exit 1
fi  

