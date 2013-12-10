#!/bin/bash


#set -x

# Parse parameters
if [ "$1" == "-t" ]
then
  shift
  TIME=$1
  shift
else
  echo "Usage: time_limit.sh -t <time> <command>"
  echo "Run <command> and kill it (and its childs) after <time> seconds."
  echo "Exit with the same exit code as the command. Code 143 if terminated by SIGTERM."
fi

# detect Cygwin
UNAME_SYSTEM=`uname -o`
if [ "${UNAME_SYSTEM}" == "Cygwin" ]
then
  HOST_OS="Cygwin"
elif [ "${UNAME_SYSTEM}" == "GNU/Linux" ]
then
  HOST_OS="Linux"
else
  HOST_OS="unkonwn"
fi  

# wait $1 seconds for process with PID $2
# return 0 if process ends before timeout $1
# return 1 if timeout is reached
function wait_for_pid {
    # wait until COMMAND is finished or timeout
    PID=$2
    TIMEOUT="$10"
    TIMER=0
    while [ ${TIMER} -lt ${TIMEOUT} ]
    do 
      if [ -e /proc/${PID} ]
      then 
        TIMER=`expr ${TIMER} + 1`
        sleep 0.1
      else
        return 0
      fi  
    done
    return 1
}



# run command on background, take its PID
#echo "$@"
"$@" &
COMMAND_PID=$!
#echo "PID: ${COMMAND_PID}"

if wait_for_pid $TIME ${COMMAND_PID}
then
  # process finished
  wait ${COMMAND_PID}
  exit $?
fi  


# kill still running COMMAND and all its childs
if [ "${HOST_OS}" == "Cygwin" ]
then
  CPIDS="${COMMAND_PID} `ps | grep "^ *[0-9]\\+ *${PID}" | awk '{print $1}'`" 
else
  CPIDS="${COMMAND_PID} $(ps --ppid ${COMMAND_PID} -o pid | tail -n +2)"
fi
for ONE_PID in ${CPIDS}
do
  kill -s SIGTERM ${ONE_PID}
  wait_for_pid 5 ${ONE_PID} &
done

# force kill for remaining
for ONE_PID in ${CPIDS}
do
  if [ -e /proc/${ONE_PID} ]
  then
    kill -9 ${ONE_PID}
  fi
done  

wait ${COMMAND_PID}
exit $?

