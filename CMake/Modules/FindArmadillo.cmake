# - Try to find Armadillo include dirs and libraries
# Usage of this module as follows:
#
# == Using Armadillo: ==
#
#   find_package( Armadillo RECQUIRED )
#   include_directories(${Armadillo_INCLUDE_DIRS})
#   add_executable(foo foo.cc)
#   target_link_libraries(foo ${Armadillo_LIBRARIES})
#
#=============================================================================
#
# This module sets the following variables:
#  Armadillo_FOUND - set to true if the library is found
#  Armadillo_INCLUDE_DIRS - list of required include directories
#  Armadillo_LIBRARIES - list of libraries to be linked 
#  Armadillo_VERSION_MAJOR - major version number
#  Armadillo_VERSION_MINOR - minor version number
#  Armadillo_VERSION_PATCH - patch version number
#  Armadillo_VERSION_STRING - version number as a string (ex: "1.0.4")
#  Armadillo_VERSION_NAME - name of the version (ex: "Antipodean Antileech")
#
#=============================================================================
# Copyright 2011 Clement Creusot
#
# This file is part of the Armadillo C++ library.
# It is provided without any warranty of fitness
# for any purpose. You can redistribute this file
# and/or modify it under the terms of the GNU
# Lesser General Public License (LGPL) as published
# by the Free Software Foundation, either version 3
# of the License or (at your option) any later version.
# (see http://www.opensource.org/licenses for more info)
#
#=============================================================================


if ( WIN32 )

  FIND_LIBRARY(Armadillo_LIBRARY
    NAMES armadillo
    PATHS "$ENV{ProgramFiles}/Armadillo/lib"  "$ENV{ProgramFiles}/Armadillo/lib64" "$ENV{ProgramFiles}/Armadillo"
    )
  FIND_PATH(Armadillo_INCLUDE_DIR
    NAMES armadillo
    PATHS "$ENV{ProgramFiles}/Armadillo/include"
    )

else ( WIN32 )  # UNIX LIKE SYSTEMS

  FIND_LIBRARY(Armadillo_LIBRARY
    NAMES armadillo
    PATHS /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib
    )
  FIND_PATH(Armadillo_INCLUDE_DIR
    NAMES armadillo
    PATHS /usr/include /usr/local/include
    )

endif ( WIN32 )


SET(Armadillo_HEADERS_FOUND FALSE)  
if(Armadillo_INCLUDE_DIR)

  # ------------------------------------------------------------------------
  #  Extract version information from <armadillo>
  # ------------------------------------------------------------------------

  # WARNING: Early releases of Armadillo didn't have the arma_version.hpp file.
  # (e.g. v.0.9.8-1 in ubuntu maverick packages (2001-03-15))
  # If the file is missing, set all values to 0  
  set(Armadillo_VERSION_MAJOR 0)
  set(Armadillo_VERSION_MINOR 0)
  set(Armadillo_VERSION_PATCH 0)
  set(Armadillo_VERSION_NAME "EARLY RELEASE")

  IF(EXISTS "${Armadillo_INCLUDE_DIR}/armadillo_bits/arma_version.hpp")

    # Read and parse armdillo version header file for version number 
    file(READ "${Armadillo_INCLUDE_DIR}/armadillo_bits/arma_version.hpp" _armadillo_HEADER_CONTENTS)
    string(REGEX REPLACE ".*#define ARMA_VERSION_MAJOR ([0-9]+).*" "\\1" Armadillo_VERSION_MAJOR "${_armadillo_HEADER_CONTENTS}")
    string(REGEX REPLACE ".*#define ARMA_VERSION_MINOR ([0-9]+).*" "\\1" Armadillo_VERSION_MINOR "${_armadillo_HEADER_CONTENTS}")
    string(REGEX REPLACE ".*#define ARMA_VERSION_PATCH ([0-9]+).*" "\\1" Armadillo_VERSION_PATCH "${_armadillo_HEADER_CONTENTS}")

    # WARNING: The number of spaces before the version name is not one.
    string(REGEX REPLACE ".*#define ARMA_VERSION_NAME\ +\"([0-9a-zA-Z\ _-]+)\".*" "\\1" Armadillo_VERSION_NAME "${_armadillo_HEADER_CONTENTS}")
  
  ENDIF(EXISTS "${Armadillo_INCLUDE_DIR}/armadillo_bits/arma_version.hpp")

  set(Armadillo_VERSION_STRING "${Armadillo_MAJOR_VERSION}.${Armadillo_MINOR_VERSION}.${Armadillo_VERSION_PATCH}")
  SET(Armadillo_INCLUDE_DIRS ${Armadillo_INCLUDE_DIR})  
  SET(Armadillo_HEADERS_FOUND TRUE)  
endif (Armadillo_INCLUDE_DIR)



#======================


IF (Armadillo_LIBRARY AND Armadillo_HEADERS_FOUND)
  SET(Armadillo_LIBRARIES ${Armadillo_LIBRARY})
  SET(Armadillo_FOUND "YES")
ELSE (Armadillo_LIBRARY AND Armadillo_HEADERS_FOUND)
  SET(Armadillo_FOUND "NO")
ENDIF (Armadillo_LIBRARY AND Armadillo_HEADERS_FOUND)

if (Armadillo_FOUND AND Armadillo_FIND_VERSION_EXACT)
  # Sets Armadillo_FOUND to TRUE if the version matchs exactly.
  set(Armadillo_FOUND FALSE)
  if(Armadillo_VERSION_MAJOR EQUAL "${Armadillo_FIND_VERSION_MAJOR}" )
    if(Armadillo_VERSION_MINOR EQUAL "${Armadillo_FIND_VERSION_MINOR}" )
      if(Armadillo_VERSION_PATCH EQUAL "${Armadillo_FIND_VERSION_PATCH}" )
	set( Armadillo_FOUND TRUE )
	set( Armadillo_EXACT_VERSION TRUE )
      endif(Armadillo_VERSION_PATCH EQUAL "${Armadillo_FIND_VERSION_PATCH}" )
    endif( Armadillo_VERSION_MINOR EQUAL "${Armadillo_FIND_VERSION_MINOR}" )
  endif( Armadillo_VERSION_MAJOR EQUAL "${Armadillo_FIND_VERSION_MAJOR}" ) 
  
  IF (Armadillo_FIND_REQUIRED)
    IF (NOT Armadillo_FOUND)
      MESSAGE(FATAL_ERROR "Could not find Armadillo Exact Version")
    ENDIF (NOT Armadillo_FOUND)
  ENDIF (Armadillo_FIND_REQUIRED)    
endif (Armadillo_FOUND AND Armadillo_FIND_VERSION_EXACT)


IF (Armadillo_FOUND)
   IF (NOT Armadillo_FIND_QUIETLY)
      MESSAGE(STATUS "Found a Armadillo library: ${Armadillo_LIBRARIES}")
   ENDIF (NOT Armadillo_FIND_QUIETLY)
ELSE (Armadillo_FOUND)
  IF (Armadillo_FIND_REQUIRED)
    IF (NOT Armadillo_HEADERS_FOUND)
      MESSAGE(FATAL_ERROR "Could not find Armadillo Headers")
    ENDIF (NOT Armadillo_HEADERS_FOUND)     
    MESSAGE(FATAL_ERROR "Could not find Armadillo library")
  ENDIF (Armadillo_FIND_REQUIRED)
ENDIF (Armadillo_FOUND)


# Hide internal variables
MARK_AS_ADVANCED(
  Armadillo_INCLUDE_DIR
  Armadillo_LIBRARY
  )



#======================