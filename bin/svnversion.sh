#!/bin/bash

#set -x

max_rev=0
for f in src/*.cc src/*.cpp include/*.h include/*.hh
do
  rev=`grep '$Revision' $f | sed 's/[^0-9]*\\([0-9]*\\).*/\\1/'` 
  if [ -z "$rev" ]; then rev=0; fi
  if [ "$rev" -gt "$max_rev" ] 
  then
    max_rev="$rev"
  fi
done
echo $max_rev
