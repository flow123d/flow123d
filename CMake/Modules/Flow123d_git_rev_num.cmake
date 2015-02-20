####
# this CMake script is called at build time to get revision of current working copy

# Arguments:
# FLOW123D_SOURCE_DIR 
#   directory of CMAKE_SOURCE_DIR
# OUTPUT_FILE_PATH
#   location where result will be saved, i.e. /path_to_file/rev_num.h
#   note:   in this action, temporary file (which is original file with suffix of .tmp 
#           i.e. OUTPUT_FILE_PATH + '.tmp') will be used 
#           after process is done, temporary file will be deleted
    

MESSAGE (STATUS "CREATING REV_NUM.H")

# get tmp file path from given FILE_PATH variable (-DFILE_PATH=path/to/file.extension)
set(OUTPUT_TMP_PATH "${OUTPUT_FILE_PATH}.tmp")

# get git-related details
include(${FLOW123D_SOURCE_DIR}/CMake/Modules/Flow123d_git_info.cmake)

# write a file with the SVN VERSION and URL define
#file(WRITE rev_num.h.tmp "#define _GIT_REVISION_ \"${GIT_DESCRIBE}\"\n#define _GIT_BRANCH_ \"${GIT_BRANCH}\"\n#define _GIT_URL_ \"${GIT_URL}\"")

configure_file(${FLOW123D_SOURCE_DIR}/CMake/rev_num_h_template ${OUTPUT_TMP_PATH})

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${OUTPUT_TMP_PATH} ${FILE_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${OUTPUT_TMP_PATH})
