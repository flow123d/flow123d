### CMake script for creating config header file
MESSAGE (STATUS "CREATING CONFIG.H")

# get tmp file path from given FILE_PATH variable (-DFILE_PATH=path/to/file.extension)
set(TMP_FILE_PATH "${FILE_PATH}.tmp")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TMP_FILE_PATH} ${FILE_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${TMP_FILE_PATH})