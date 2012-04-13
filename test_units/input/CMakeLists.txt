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


if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-std=gnu++0x")        ## Use C++11 for auto
    #define_test(enum_macro)
    
    #define_test(input_interface)
    
    define_test(input_type)

    #set(CMAKE_CXX_FLAGS "-std=gnu++0x")        ## Use C++11 for multiline raw strings
    #define_test(generic_node)
    
endif()





