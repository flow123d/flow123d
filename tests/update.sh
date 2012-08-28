#!/bin/bash

# Syntax:
#  update.sh  dir_to_update  source_dir
#
#  Recursively updates all files in <dir_to_update> by corresponding files in <source_dir>
#  Reports an error if the file or irectory in <source_dir> doesn't exist.

#set -x

function update_dir {
  
  local target_dir=$1
  local source_dir=$2
  
  for item in `ls ${target_dir}/`
  do
    if [ -d "${target_dir}/${item}" ]
    then 
      if [ ! -d "${source_dir}/${item}" ]
      then
        echo "Missing directory in test result."
        return 1
      fi  
      update_dir "${target_dir}/${item}" "${source_dir}/${item}"
    else
      if [ ! -e "${source_dir}/${item}" ]
      then
        echo "Missing file '${source_dir}/${item}' in test result."
        return 1
      fi
      
      echo "Copping ${source_dir}/${item} -> ${target_dir}/${item}"
      cp ${source_dir}/${item} ${target_dir}/${item}
    fi
  done
}

update_dir $1 $2
