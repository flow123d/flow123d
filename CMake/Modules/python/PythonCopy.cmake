# PythonCopy module
#  - this module will copy Python interpret and its core packages into
#    PYTHON_COPY_DESTINATION folder
#    
#    
#  - Following variables MUST BE SET:
#    - PYTHON_ROOT
#           - absolute path to python root e.g. /usr/lib/python2.7
#           
#    - PYTHON_COPY_PACKAGES 
#           - list of all core modules which will be copied
#             each entry represents folder name realtive to root
#             e.g. set(PYTHON_COPY_MODULES encodings json)
#    - PY_BUILD_PREFIX
#           - directory where libs will be copied
#             in ${PY_BUILD_PREFIX}/lib/python2.7
#             
#   
#   - this module export following variables:
#           - PYTHON_COPY_DESTINATION
#               - location where all python files are copied
#           - PYTHON_COPY_INCLUDE_TARGET
#               - location where python h files for embedded python will be copied
# ------------------------------------------------------------------------------


# check root presence
if(NOT EXISTS ${PYTHON_ROOT})
    message(FATAL_ERROR "Invalid PYTHON_ROOT. Folder ${PYTHON_ROOT} does not exists.")
endif()
    
# set extra module path to binary lib dir
set(PYTHON_COPY_DESTINATION ${PY_BUILD_PREFIX}/lib/python2.7)

# get files and folders from root (non recursive)
file(GLOB TMP ${PYTHON_ROOT}/*)
# go through all files and copy them to root
foreach(ITEM ${TMP})
    if(NOT IS_DIRECTORY ${ITEM})
        if(EXISTS ${ITEM})
            file(COPY ${ITEM} DESTINATION ${PYTHON_COPY_DESTINATION})
        endif()
    endif()
endforeach()

# process all core modules
foreach(ITEM ${PYTHON_COPY_PACKAGES})
    if(EXISTS ${PYTHON_ROOT}/${ITEM})
        file(COPY ${PYTHON_ROOT}/${ITEM} DESTINATION ${PYTHON_COPY_DESTINATION})
    endif()
endforeach()

# copy flow123d python scripts and modules
execute_process(COMMAND
    cp -rT ${CMAKE_SOURCE_DIR}/src/python ${PY_BUILD_PREFIX}/lib/flow123d)

# copy python include (h files for embedded python)
set(PYTHON_COPY_INCLUDE_TARGET include/)
execute_process(COMMAND
    mkdir -p ${PY_BUILD_PREFIX}/include)
# recursive copy, not overwrite existing files
# we copy all and exclude some paths and files during install    
execute_process(COMMAND
    cp -rf ${PYTHON_INCLUDE_DIR} ${PY_BUILD_PREFIX}/include)