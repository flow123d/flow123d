####
# this CMake script is called at build time to get revision of current working copy

# the FindSubversion.cmake module is part of the standard distribution
include(FindSubversion)

if (Subversion_FOUND) 
  # extract working copy information for SOURCE_DIR into Flow variable
  Subversion_WC_INFO(${SOURCE_DIR} Flow)
  message(STATUS "SVN working copy revision: ${Flow_WC_REVISION}")
  message(STATUS "SVN working copy URL: ${Flow_WC_URL}")
else ()
  set(Flow_WC_REVISION "(unknown revision)")
endif()

# write a file with the SVN VERSION and URL define
file(WRITE rev_num.h.tmp "#define _PROGRAM_REVISION_ \"${Flow_WC_REVISION}\"\n#define _PROGRAM_BRANCH_ \"${Flow_WC_URL}\"\n")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different rev_num.h.tmp rev_num.h)
#execute_process(COMMAND ${CMAKE_COMMAND} -E remove rev_num.h.tmp)