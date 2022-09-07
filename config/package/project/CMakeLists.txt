# project name
project(Flow123d C)

# minimum version
cmake_minimum_required(VERSION 2.8.8)


# configurable via makefile
SET(IMAGE_TAG           "flow123d/v200" CACHE STRING "Docker image name (how it will be tagged)")
SET(DEST                "destination/folder" CACHE STRING "Destination folder, where tests and doc and images are stored")
SET(FLOW123D_ROOT       "/opt/flow123d/flow123d" CACHE STRING "Location of the Flow123d repository root")
SET(FLOW_VERSION        "1.0.0" CACHE STRING "Dot separated")


message(STATUS "Packaging Flow123d version ${FLOW_VERSION}")
message(STATUS "Flow123d root: '${FLOW123D_ROOT}'")


  message(STATUS "Creating Linux package")

  # configure script so that IMAGE_TAG is correctly set
  configure_file(src/linux/install.sh     pack/install.sh     @ONLY)
  configure_file(src/linux/uninstall.sh   pack/uninstall.sh   @ONLY)
  configure_file(src/linux/flow123d.sh    pack/flow123d.sh    @ONLY)
  configure_file(src/linux/fterm.sh       pack/fterm.sh       @ONLY)
  configure_file(src/linux/runtest.sh     pack/runtest.sh     @ONLY)

 
  # install tests and docs
  INSTALL(DIRECTORY ${DEST}/htmldoc/                          DESTINATION doc/htmldoc/)
  INSTALL(FILES     ${DEST}/flow123d_${FLOW_VERSION}_doc.pdf  DESTINATION doc/)

  INSTALL(DIRECTORY ${DEST}/tests/                        DESTINATION tests/)
  INSTALL(PROGRAMS  ${CMAKE_BINARY_DIR}/pack/runtest.sh   DESTINATION tests/)

  # install and uninstall scripts to root
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/pack/install.sh    DESTINATION ./)
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/pack/uninstall.sh  DESTINATION ./)
  # add shortcuts to bin
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/pack/fterm.sh      DESTINATION bin/)
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/pack/flow123d.sh   DESTINATION bin/)
  INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/pack/runtest.sh    DESTINATION bin/)

  # install configure
  INSTALL(PROGRAMS ${DEST}/bin/fterm                      DESTINATION bin
                                                          RENAME      configure)

  SET(CPACK_PACKAGE_FILE_NAME "flow123d_${FLOW_VERSION}")
  SET(CPACK_GENERATOR "TGZ")

# create package rule
include(CPack)
