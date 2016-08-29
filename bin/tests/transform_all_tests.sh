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

  if [ -f ${f%.yaml}.con ]  # skip *.yaml converted from *.con
  then
    continue
  fi
  if [ ! "${f%.orig.yaml}" == "$f" ]   # skip *.orig.yaml  
  then
    continue
  fi

  f_base=${f%.yaml}  
  for forig in ${f_base}.*.orig.yaml  # revert last transform of yaml before new transform
  do
    if [ "${forig}" == "${f_base}.0.orig.yaml" ]
    then
        cp -f ${forig} ${f}
    fi
    rm ${forig}    
  done

  #echo $f
  python3 "${0%/*}/../input_convert.py" $f
done
