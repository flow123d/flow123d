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
echo "$@"
$@ &
PID=$!

(
  sleep $TIME; 
  # kill only if the program still runs
  ps cax | grep "^ *${PID}" > /dev/null;
  if [ $? -eq 0 ]; 
  then  kill ${PID};
  fi
)&

# wait until the program finish
wait $PID
