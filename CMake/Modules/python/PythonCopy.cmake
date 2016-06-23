# PythonCopy module
#  - this module will copy Python interpret and its core packages into
#    PYTHON_COPY_DESTINATION folder
#    
#    
#  - Following variables MUST BE SET:
#    - PYTHON_COPY_ROOT
#           - absolute path to python root e.g. /usr/lib/python2.7
#           
#    - PYTHON_COPY_PACKAGES 
#           - list of all core modules which will be copied
#             each entry represents folder name realtive to root
#             e.g. set(PYTHON_COPY_MODULES encodings json)
#             
#   - PYTHON_COPY_DIST_PACKAGES
#           - list of all 3rd party packages (installed using pip or other tool)
#             each entry is either
#                - relative path:
#                   - relativity is taken towards PYTHON_COPY_ROOT
#                   - all relativity is kept during copy so dist-packages/yaml
#                     will create folder dist-packages/yaml
#                - absolute path:
#                   - only last 2 folders are taken into account for relativity
#                     (we cannot conclusively estimate where relativity begins)
#                     e.g. /usr/local/lib/dist-packages/markdown
#                     will create folder dist-packages/markdown
#   
#   
#   - this module export following variables:
#           - PYTHON_COPY_DESTINATION
#               - location where all python files are copied
#           - PYTHON_COPY_INCLUDE_TARGET
#               - location where python h files for embedded python will be copied
# ------------------------------------------------------------------------------


if(NOT DEFINED PYTHON_COPY_ROOT)
    message(STATUS "--------------------------------------------------")
    message(STATUS "PYTHON_COPY_ROOT not set! Using default configuration which may not work on this system")
    message(STATUS "--------------------------------------------------")
    
    set(PYTHON_COPY_ROOT /usr/lib/python2.7)
    set(PYTHON_COPY_PACKAGES
        bsddb compiler ctypes curses distutils email encodings hotshot importlib 
        json lib-dynload lib-tk lib2to3 logging multiprocessing 
        plat-x86_64-linux-gnu pydoc_data sqlite3 test unittest wsgiref xml)
    set(PYTHON_SYSPATH
        plat-x86_64-linux-gnu lib-tk lib-old lib-dynload dist-packages site-packages)
    set(PYTHON_COPY_DIST_PACKAGES)
endif()

# check root presence
if(NOT EXISTS ${PYTHON_COPY_ROOT})
    message(FATAL_ERROR "Invalid PYTHON_COPY_ROOT. Folder ${PYTHON_COPY_ROOT} does not exists.")
endif()
    
# set extra module path to binary lib dir
set(PYTHON_COPY_DESTINATION ${PY_BUILD_PREFIX}/lib/python2.7)

# get files and folders from root (non recursive)
file(GLOB TMP ${PYTHON_COPY_ROOT}/*)

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
    if(EXISTS ${PYTHON_COPY_ROOT}/${ITEM})
        file(COPY ${PYTHON_COPY_ROOT}/${ITEM} DESTINATION ${PYTHON_COPY_DESTINATION})
    endif()
endforeach()

# process dist-packages
foreach(ITEM ${PYTHON_COPY_DIST_PACKAGES})
    if(NOT IS_ABSOLUTE ${ITEM})
        set(ITEM_ABSPATH ${PYTHON_COPY_ROOT}/${ITEM})
        set(ITEM_PATH ${ITEM})
    else()
        set(ITEM_ABSPATH ${ITEM})
        string(REGEX REPLACE "^.*/(.+/.+)$" "\\1" ITEM_PATH ${ITEM_ABSPATH})
    endif()
    
    string(REGEX REPLACE "(.+)/(.+)$" "\\1" PACKAGE_PATH ${ITEM_PATH})
    if(EXISTS ${ITEM_ABSPATH})
        file(COPY ${ITEM_ABSPATH} DESTINATION ${PYTHON_COPY_DESTINATION}/${PACKAGE_PATH})
    endif()
endforeach()

# copy flow123d python scripts and modules
execute_process(COMMAND
    cp -rT ${CMAKE_SOURCE_DIR}/src/python ${PY_BUILD_PREFIX}/lib/flow123d)

# copy python include (h files for embedded python)
# TODO: where is PYTHON_SUBDIR variable set?
set(PYTHON_COPY_INCLUDE_TARGET include/${PYTHON_SUBDIR})
execute_process(COMMAND
    mkdir -p ${PY_BUILD_PREFIX}/include)
# recursive copy, not overwrite existing files
# we copy all and exclude some paths and files during install    
execute_process(COMMAND
    cp -rn ${PYTHON_INCLUDE_DIR} ${PY_BUILD_PREFIX}/include)