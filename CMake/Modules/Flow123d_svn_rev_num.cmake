####
# this CMake script is called at build time to get revision of current working copy

# the FindSubversion.cmake module is part of the standard distribution
include(FindSubversion)

if (Subversion_FOUND) 
  # extract working copy information for SOURCE_DIR into Flow variable
  Subversion_WC_INFO(${SOURCE_DIR} Flow)
  message(STATUS "SVN revision: ${Flow_WC_REVISION}")
else ()
  set(Flow_WC_REVISION "unknown revision")
endif()

# write a file with the SVNVERSION define
file(WRITE rev_num.h.tmp "#define REVISION \"${Flow_WC_REVISION}\"\n")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different rev_num.h.tmp rev_num.h)
#execute_process(COMMAND ${CMAKE_COMMAND} -E remove rev_num.h.tmp)