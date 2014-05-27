#!/bin/bash

cd ${0%/*}/../build_tree/bin
if [ -x  ./flow123d.exe ]
then
  libs=`cygcheck.exe ./flow123d.exe | grep ".*\cyg.*.dll"`
  echo "Copy flow123d.ex libraries ..."
  cp $libs .
else
  echo "No binary. No libraries."
fi  
