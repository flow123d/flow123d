#!/bin/bash

# This script is called during installation phase.
# Copy shared libraries (on Linux) or DLLs (under Cygwin) on which 
# flow123d binary depends in order to have self-containd package.
#
# Shared libraries are placed to ../lib relative to executable, while
# the same path is set to RPATH of the executable (CMake do it just during installation phase, 
# the executable is linked with full paths and relinked during installation)
# We won't to have provided libraries in the directory of executable since it usually breaks
# function of system command when called from the directory that contains the libraries.
#
# DLLs are placed in the same directory executable. (There is no RPATH under windows.)

set -x


cd ${0%/*}/../build_tree/bin
if [ "$1" == "cygwin" ]
then
  if [ -x  ./flow123d.exe ]
  then
    libs=`cygcheck.exe ./flow123d.exe | grep ".*\cyg.*.dll"`

    echo "Copy flow123d shared libraries ..."
    echo $libs
    cp $libs .

  fi
else  
  if [ -x  ./flow123d ]
  then
    # detect shared libs
    loader=`ldd flow123d | grep "/ld" | sed 's/^\t*\(\/.*\) (.*)$/\1/'`
    libs=`ldd flow123d | grep -n "^/" | sed 's/^.* => \(\/.*\) (.*)$/\1/'`
    
    # copy shared libs
    echo "Copy flow123d shared libraries ..."
    if [ ! -d ../lib ]; then mkdir ../lib; fi
    cp ${libs} ../lib
    
    # manually copy and rename file to bin dir so cmake can make it executable
    cp ${loader} ./ld-linux-loader.so
    chmod +x "./ld-linux-loader.so"
  fi
fi

#if [ -n "${libs}" ]
#then
#else
#  echo "No binary. No libraries."
#fi  
