####
# this CMake script is called at build time to get revision of current working copy

include(FindGit)

if (GIT_FOUND) 
  
  message(STATUS "PWD: ${SOURCE_DIR}")
  # get GIT branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Get the describe string - tag + num. of commits since tag + short hash
  execute_process(
    COMMAND git describe
    WORKING_DIRECTORY ${SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DESCRIBE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  
  # Get remote branch we are tracking
  execute_process(
    COMMAND git branch --list -vv ${GIT_BRANCH}
    WORKING_DIRECTORY ${SOURCE_DIR}
    OUTPUT_VARIABLE GIT_FULL_BRANCH_INFO
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  STRING(REGEX REPLACE ".* \\[([^]]*)\\] .*" "\\1" GIT_REMOTE_BRANCH "${GIT_FULL_BRANCH_INFO}")
  STRING(REGEX REPLACE "([^\\w]*)/.*" "\\1" GIT_REMOTE "${GIT_REMOTE_BRANCH}")
  
  # Get Fetch URL of the remote
  execute_process(
    COMMAND git remote show -n ${GIT_REMOTE}
    WORKING_DIRECTORY ${SOURCE_DIR}
    OUTPUT_VARIABLE GIT_FULL_URL_INFO
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  STRING(REGEX REPLACE "\"" "\\\\\"" GIT_FULL_URL_ESC_QUOTES "${GIT_FULL_URL_INFO}")
  STRING(REGEX MATCH "Fetch URL:[^\n]*" GIT_URL_TMP "${GIT_FULL_URL_ESC_QUOTES}")
  STRING(REGEX REPLACE "Fetch URL:([^\n]*)" "\1" GIT_URL "${GIT_URL_TMP}")
  
  
  message(STATUS "GIT_REMOTE: ${GIT_REMOTE}")
  message(STATUS "GIT_BRANCH: ${GIT_BRANCH}")
  message(STATUS "GIT_DESCRIBE: ${GIT_DESCRIBE}")
  message(STATUS "GIT_REMOTE_BRANCH: ${GIT_REMOTE_BRANCH}")
  message(STATUS "GIT_FULL_URL: ${GIT_FULL_URL_INFO}")
  message(STATUS "GIT_URL: ${GIT_URL}")
endif()

# write a file with the SVN VERSION and URL define
file(WRITE rev_num.h.tmp "#define _GIT_REVISION_ \"${GIT_DESCRIBE}\"\n#define _GIT_BRANCH_ \"${GIT_BRANCH}\"\n#define _GIT_URL_ \"${GIT_URL}\"")

# copy the file to the final header only if the version changes
# reduces needless rebuilds
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different rev_num.h.tmp rev_num.h)
#execute_process(COMMAND ${CMAKE_COMMAND} -E remove rev_num.h.tmp)

set(GIT_DESCRIBE ${GIT_DESCRIBE} CACHE INTERNAL "Human readable description of last git commit.")