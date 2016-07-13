#!/bin/bash
#
# Convert all CON files in directory "tests" to YAML format.
# Also apply possible transformations due to the change in input source tree.
# The transformations have to be defined in the file
# "src/python/GeoMop/ModelEditor/resources/transformations/${TRANSFORM_FILE}.json".


# Relative path to "tests" directory from directory,
# where this script is placed
TESTS_DIR="../../tests"
# Relative path to "tests" directory from current/working directory
TESTS_DIR="${0%/*}/${TESTS_DIR}"


for f in $TESTS_DIR/*/*.con $TESTS_DIR/*/*.yaml
do

  if [ -f ${f%.yaml}.con ]
  then
    continue
  fi
  if [ ! "${f%.orig.yaml}" == "$f" ]
  then
    continue
  fi

  #echo $f
  "${0%/*}/transform_test.sh" $f
done
