####
# this CMake script is called at build time to get revision of current working copy

# the FindSubversion.cmake module is part of the standard distribution
include(FindGit)

if (GIT_FOUND) 
  
  # get GIT branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the describe string - tag + num. of commits since tag + short hash
  execute_process(
    COMMAND git describe
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  
  
  message(STATUS "GIT_BRANCH: ${GIT_BRANCH}")
  message(STATUS "GIT_DESCRIBE: ${GIT_DESCRIBE}")
endif()

# write a file with the SVN VERSION and URL define
file(WRITE rev_num.h.tmp "#define _PROGRAM_REVISION_ \"${GIT_DESCRIBE}\"\n#define _PROGRAM_BRANCH_ \"${GIT_BRANCH}\"\n")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different rev_num.h.tmp rev_num.h)
#execute_process(COMMAND ${CMAKE_COMMAND} -E remove rev_num.h.tmp)