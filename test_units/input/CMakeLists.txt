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
# $Id: CMakeLists.txt 1567 2012-02-28 13:24:58Z jan.brezina $
# $Revision: 1567 $
# $LastChangedBy: jan.brezina $
# $LastChangedDate: 2012-02-28 14:24:58 +0100 (Tue, 28 Feb 2012) $
#

#set(CMAKE_INCLUDE_CURRENT_DIR ON)



add_custom_target( input_test )
set(super_test "input_test")
set(super_test_libs input_lib)

# for given class_name assume class_name_test.cpp source
# and class_name_test binary
# and class_name_test test target
# use libraries in ${super_test_libs}
#
# I want to run tests from source dir, but since there is no portable way 
# how to make link from source dir to Makefile generated in bulid dir 
# I let main makefile to make such links.
macro(define_test class_name)
  set(test_source "${class_name}_test.cpp")
  set(test_binary "${class_name}_test_bin")
  set(test_name "${class_name}_test")
  
  # executable rule
  add_executable(${test_binary} EXCLUDE_FROM_ALL  ${test_source}) 
  target_link_libraries(${test_binary} ${super_test_libs} gtest gtest_main)

  # rule to run the test
  add_custom_target(${test_name} 
      COMMAND ${test_binary}
      DEPENDS ${test_binary}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # register the test into super_test
  add_dependencies(${super_test} ${test_name})

  # let ctest call make with appropriate  target
  add_test(NAME ${test_name} COMMAND ${CMAKE_MAKE_PROGRAM} ${test_name})
#  add_test(NAME ${test_binary} COMMAND 
#    ${CMAKE_CTEST_COMMAND}
#    --build-and-test "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}"
#    --build-nocmake
#    --build-generator ${CMAKE_GENERATOR}
#    --build-makeprogram ${CMAKE_MAKE_PROGRAM}
#    --build-target ${test_binary}
#    --test-command ${test_binary}
#    --verbose)
endmacro()

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=gnu++0x")        ## Optimize
    define_test(generic_node)
endif()





