#!/usr/bin/env bash
#
# This script downloads and runs flow123d

# check if stdout is a terminal...
if [[ -z "$nocolor" ]]; then
  if test -t 1; then
      # see if it supports colors...
      ncolors=$(tput colors)
      if test -n "$ncolors" && test $ncolors -ge 8; then
          bold="$(tput bold)"
          reset="$(tput sgr0)"
          red="$(tput setaf 1)"
          green="$(tput setaf 2)"
          yellow="$(tput setaf 3)"
          blue="$(tput setaf 4)"
          bblue="$bold$blue"
          bgreen="$bold$green"
          byellow="$bold$yellow"
          bred="$bold$red"
      fi
  fi
fi

function dbg() {
  if [[ $verbose -eq 1 ]]; then
    echo -e "$bgreen[DBG]$reset $@"
  fi
}
function dbg2() {
  if [[ $verbose -eq 1 ]]; then
    echo -e "$byellow[DBG]$reset $@"
  fi
}
function dbgc() {
  if [[ $verbose -eq 1 ]]; then
    echo -e "$bblue[RUN]$reset $@"
  fi
  $@
}

version=@IMAGE_TAG@
if [[ -z "$version" ]]; then
    echo -e "$byellow no docker image version specified$reset"
    exit 1
fi

mnt=$(pwd)
verbose=1
privileged=0
ACTION=shell
cwd=$(pwd)
uid=$(id -u)
gid=$(id -g)
uname=flow
tty="true"
interactive="true"

# define usage
function usage() {
  cat << EOF
usage: flow123d [--help] [--privileged] [--version <VERSION>] [--mount <MOUNT>] [--workdir <WORKDIR>] [<ACTION>] [<args>]
Where ACTION can be:
  shell                 (default behaviour) Enter interactive shell
  run                   Execute flow123d and pass all given arguments
  help                  Print this message and exits

  -v, --version <VERSION>
      Version which will be used, default value is ${bblue}$version${reset}

  -m, --mount <MOUNT>
      A directory which will be mounted (files thus will be accessible),
        default value is a current directory ${bgreen}\`pwd\`${reset}
        current value              ${bblue}$mnt${reset}

  -w, --workdir <WORKDIR>
      A working directory, default value is current working directory,
        current value              ${bblue}$cwd${reset}

  -p, --privileged ${reset}
      Will add --privileged=true when starting docker container,
      this options carries a security risk but should deal with SELinux mounting
      issues

  --tty true|false ${reset}
      Will set tty mode if true non-tty or false
      current value                ${bblue}--tty $tty${reset}

  --interactive true|false ${reset}
      Will set interactive mode if true non-interactive or false
      current value                ${bblue}--interactive $interactive${reset}

  <args>
      Additional arguments which are passed to the flow123d (in case ACTION is run)
      otherwise passed to the docker run

  -h, --help
      Print this message and exits
EOF
}

while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
      -h|--help|help)
        usage
        exit 0
      ;;
      -m|--mount)
        mnt="$2"
        shift # past argument
        shift # past value
      ;;
      -w|--workdir)
        cwd="$2"
        shift # past argument
        shift # past value
      ;;
      -v|--version)
        version="$2"
        shift # past argument
        shift # past value
      ;;
      -p|--privileged)
        privileged=1
        shift
      ;;
      --tty)
        tty="$2"
        shift
        shift
      ;;
     --interactive)
        interactive="$2"
        shift
        shift
      ;;
      run|shell|exec)
        ACTION="$1"
        shift # past argument
        
        # in case someone did the fterm run -- <args>
        # we strip the first double dashes
        if [[ "$1" == "--" ]]; then
          shift
        fi
        break
      ;;
      *)
        break
      ;;
  esac
done

# remove flow123d prefix
final_version=${version#flow123d/}


if [[ $privileged == "1" ]]; then
  priv_true="--privileged=true"
fi


flags="-e uid=$uid -e gid=$gid -ewho=$uname -v $mnt:$mnt -w $cwd $priv_true"


if [[ "$ACTION" == "shell" ]]; then
  dbgc docker run --tty=$tty --interactive=$interactive --rm $flags flow123d/$final_version

elif [[ "$ACTION" == "run" ]]; then
  dbgc docker run --tty=$tty --interactive=$interactive --rm $flags flow123d/$final_version flow123d "$@"

elif [[ "$ACTION" == "exec" ]]; then
  dbgc docker run --tty=$tty --interactive=$interactive --rm $flags flow123d/$final_version "$@"
fi

exit $?
