#!/bin/bash
#
# Usage: transform_test.sh file1.con file2.con ...
#
# Convert CON file(s) to YAML format.
# Also apply possible transformations due to the change in input source tree.
# The transformations have to be defined in the file
# "src/python/GeoMop/ModelEditor/resources/transformations/main.json".

# Relative path to "importer.py" script from directory,
# where this script is placed
IMPORTER_PY="../../src/python/GeoMop/ModelEditor/importer.py"
# Relative path to "importer.py" script from current/working directory
IMPORTER_PY="${0%/*}/${IMPORTER_PY}"


for f in $@; do
  echo "Processing $f";
  `rm -f ${f%.con}.yaml`;
  `python3 $IMPORTER_PY --transformation-name main --con_file $f`;
done