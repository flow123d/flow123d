#!/bin/bash
#
# This script is used in Jenkins publish job to collect various files we want 
# to publish into common directory.  This directory is named according to the build version
# (given by GIT_VERSION_FULL cmake variable extracted from CMakeCache.txt). Furthermore, we rename packages
# to canonical form: flow123d_<VERSION>_<PLATFORM>.<SUFFIX>.
#
# Script do not copy artifacts that are not present.

set -x

BIN_DIR="${0%/*}"
BUILD_TREE="${BIN_DIR}/../build_tree" 
SOURCE_DIR="${BIN_DIR}/.."

function get_cmake_var {
  VAR_NAME=$1
  VALUE=`grep "$1" "${BUILD_TREE}/CMakeCache.txt" | sed 's/^.*=//'`
  eval $VAR_NAME=\$VALUE
}

get_cmake_var GIT_VERSION_FULL
get_cmake_var PLATFORM_NAME

PUBLISH_DIR="${BUILD_TREE}/publish_dir/${GIT_VERSION_FULL}"
FLOW_BIN="${BUILD_TREE}/bin/flow123d"

function safe_copy {
  TARGET=${PUBLISH_DIR}
  if [ -n "$2" ]
  then 
    TARGET="$2"
  fi  
  
  if [ -e "$1" ]
  then
    cp -r "$1" "${TARGET}"
  fi  
}

rm -rf "${PUBLISH_DIR}"
mkdir -p "${PUBLISH_DIR}"


#${BIN_DIR}/md_convert.sh ${BIN_DIR}/../README.md
#${BIN_DIR}/md_convert.sh ${BIN_DIR}/../CHANGES.md
#safe_copy "${BUILD_TREE}/README.html"
#safe_copy "${BUILD_TREE}/CHANGES.html"
safe_copy "${SOURCE_DIR}/CHANGES.md"
safe_copy "${SOURCE_DIR}/README.md"
safe_copy "${BUILD_TREE}/doc/reference_manual/flow123d_doc.pdf"
safe_copy "${BUILD_TREE}/doc/online-doc/flow123d" "${PUBLISH_DIR}/source_doc"
safe_copy "${BUILD_TREE}/Flow123d-${GIT_VERSION_FULL}-Linux.tar.gz" \
          "${PUBLISH_DIR}/flow123d_${GIT_VERSION_FULL}_${PLATFORM_NAME}.tar.gz"
safe_copy "${BUILD_TREE}/Flow123d-${GIT_VERSION_FULL}-Source.tar.gz" \
          "${PUBLISH_DIR}/flow123d_${GIT_VERSION_FULL}_source.tar.gz"
#safe_copy "${BUILD_TREE}/Flow123d-${GIT_VERSION_FULL}-CYGWIN.zip" \
#          "${PUBLISH_DIR}/flow123d_${GIT_VERSION_FULL}_${PLATFORM_NAME}.zip"
safe_copy "${BUILD_TREE}/_CPack_Packages/CYGWIN/NSIS/Flow123d-${GIT_VERSION_FULL}-CYGWIN.exe" \
          "${PUBLISH_DIR}/flow123d_${GIT_VERSION_FULL}_${PLATFORM_NAME}.exe"

