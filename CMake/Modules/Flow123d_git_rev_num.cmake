####
# this CMake script is called at build time to get revision of current working copy

include(${FLOW123D_SOURCE_DIR}/CMake/Modules/Flow123d_git_info.cmake)

# write a file with the SVN VERSION and URL define
#file(WRITE rev_num.h.tmp "#define _GIT_REVISION_ \"${GIT_DESCRIBE}\"\n#define _GIT_BRANCH_ \"${GIT_BRANCH}\"\n#define _GIT_URL_ \"${GIT_URL}\"")

configure_file(${FLOW123D_SOURCE_DIR}/CMake/rev_num_h_template ${CMAKE_BINARY_DIR}/rev_num.h.tmp)

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_BINARY_DIR}/rev_num.h.tmp rev_num.h)
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/rev_num.h.tmp)
