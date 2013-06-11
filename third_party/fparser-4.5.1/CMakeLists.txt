# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Centre: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License version 3 as published by the Free Software Foundation.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program; if not,
# write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 021110-1307, USA.
#
# $Id: CMakeLists.txt 1306 2011-08-30 13:53:56Z jiri.jenicek $
# $Revision: 1306 $
# $LastChangedBy: jiri.jenicek $
# $LastChangedDate: 2011-08-30 15:53:56 +0200 (Tue, 30 Aug 2011) $
#

message(STATUS "CMakeLists.txt - FParser")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# select source files that compose whole flow123d library
# todo: split into parts, e.g. system, output, input, mesh, darcy_mh, transport, .. 
# every library then can be testetd independently
set(sources_lib_fparser
  fpoptimizer.cc
  fparser.cc
)

# add tager for libfparser
add_library(fparser ${sources_lib_fparser})
set_target_properties(fparser PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib
)



