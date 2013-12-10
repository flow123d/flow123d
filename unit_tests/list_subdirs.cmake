#  List subdirectories
#  -------------------
#  accepted parameters:
#   DIRECTORY - directory to list  
#   OUTPUT - file to which write the list (overwritten)
# CMake script that list subdirectories (relative) in the DIRECTORY and writes the result too stdout.

   file(REMOVE ${OUTPUT}) 
   
   file(GLOB dir_list ${DIRECTORY}/* )
   list(SORT dir_list)
   #message("LIST: ${dir_list}")
   set(list_of_dirs "")
   foreach(list_item ${dir_list})
     if(IS_DIRECTORY ${dir}/${list_item})
         string(REGEX REPLACE ".*/([^/]*)" "\\1"
                sub_dir 
                ${list_item}
                )
         file(APPEND ${OUTPUT} "${sub_dir}\n")       
     endif()
   endforeach()
