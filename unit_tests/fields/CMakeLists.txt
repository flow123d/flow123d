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

set(libs coupling_lib mesh_lib fem_lib flow123d_system ${Armadillo_LIBRARIES} ${Armadillo_LINK_LIBRARIES} ${PYTHON_LIBRARIES})
add_test_directory("${libs}")



define_test(field_value)

define_mpi_test(field_algo_base 1)
define_mpi_test(field_algo_base 2)

define_mpi_test(field_set 1)

define_mpi_test(field_const 1)
define_mpi_test(field_formula 1)
define_mpi_test(field_python 1)
define_mpi_test(field_fe 1)
define_mpi_test(multi_field 1)
define_mpi_test(field_time_function 1)

define_mpi_test(generic_field 1)

#define_mpi_test(field_speed 1 30)
#define_mpi_test(field_constant_speed 1)
#define_mpi_test(field_fe_speed 1)
#define_mpi_test(field_model_speed 1)
       
define_mpi_test(eval_subset 1)
define_mpi_test(field_value_cache 1)

define_test(bparser)
define_test(field_model)




