#!/bin/bash
#
# Usage: transform_test.sh file1.con file2.con ...
#
# Convert CON file(s) to YAML format.
# Also apply possible transformations due to the change in input source tree.
# The transformations have to be defined in the file
# "src/python/GeoMop/ModelEditor/resources/transformations/${TRANSFORM_FILE}.json".

TRANSFORM_FILE_CON=flow123d_1.8.6_to_2.0.0
TRANSFORM_FILE_YAML=flow123d_2.0.0_rc_to_2.0.0

# Relative path to "importer.py" script from directory,
# where this script is placed
IMPORTER_PY="../yaml_converter/yaml_converter.py"
# Relative path to "importer.py" script from current/working directory
IMPORTER_PY="${0%/*}/${IMPORTER_PY}"


for f in $@ 
do
  echo "Processing $f"
  if [ "${f%.yaml}" == "$f" ]
  then 
        # con file
        echo "Conversion of CON files not available." 
        # rm -f "${f%.con}.yaml"
        # python3 "$IMPORTER_PY" --transformation-name "${TRANSFORM_FILE_CON}" --con_file "$f"
  else     
        # yaml file
        
        new_name="${f%.yaml}.orig.yaml"
        if [ ! -f "${new_name}" ]
        then
            cp "$f" "${new_name}"        
        fi    
        rm "$f"
        python3 "$IMPORTER_PY"  "${new_name}"
        converted=${new_name%.yaml}.new.yaml
        if [ -f $converted ]
            mv ${converted} $f
        else
            mv $new_name $f
        fi    
  fi  
done
