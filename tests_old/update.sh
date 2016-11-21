#!/bin/bash

# Syntax:
#  update.sh  source_dir dir_to_update  
#
#  Recursively updates all files in <dir_to_update> by corresponding files in <source_dir>
#  Reports an error if the file or irectory in <source_dir> doesn't exist.

#set -x

function update_dir {
  
  local target_dir=$2
  local source_dir=$1
  
  for item in `ls ${target_dir}/`
  do
    if [ -d "${target_dir}/${item}" ]
    then 
      dir_item=
      if [ -d "${source_dir}/${item}" ]
      then
        update_dir "${source_dir}/${item}" "${target_dir}/${item}"
      elif [ -d "${source_dir}/${item}.1" ]  
      then
        update_dir "${source_dir}/${item}.1" "${target_dir}/${item}"
      else  
        echo "Missing source directory '${source_dir}/${item}.*' "
        return 1
      fi  
      
    else
      if [ ! -e "${source_dir}/${item}" ]
      then
        echo "Missing source file '${source_dir}/${item}' "
        return 1
      fi
      
      echo "Copying ${source_dir}/${item} -> ${target_dir}/${item}"
      cp ${source_dir}/${item} ${target_dir}/${item}
    fi
  done
}

update_dir $1 $2
