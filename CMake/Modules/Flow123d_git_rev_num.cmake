####
# this CMake script is called at build time to get revision of current working copy
MESSAGE (STATUS "CREATING REV_NUM.H")

# get tmp file path from given FILE_PATH variable (-DFILE_PATH=path/to/file.extension)
set(TMP_FILE_PATH "${FILE_PATH}.tmp")

# get git-related details
include(${FLOW123D_SOURCE_DIR}/CMake/Modules/Flow123d_git_info.cmake)

# write a file with the SVN VERSION and URL define
#file(WRITE rev_num.h.tmp "#define _GIT_REVISION_ \"${GIT_DESCRIBE}\"\n#define _GIT_BRANCH_ \"${GIT_BRANCH}\"\n#define _GIT_URL_ \"${GIT_URL}\"")

configure_file(${FLOW123D_SOURCE_DIR}/CMake/rev_num_h_template ${TMP_FILE_PATH})

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${TMP_FILE_PATH} ${FILE_PATH})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${TMP_FILE_PATH})
