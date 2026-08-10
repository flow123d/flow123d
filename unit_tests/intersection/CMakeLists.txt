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

#set(super_test_libs fast_intersection_lib field_lib mesh_lib system_lib)

set(libs mesh_lib flow123d_system io_lib)
add_test_directory("${libs}")
    
    define_mpi_test(compute_intersection_12 1)
    define_mpi_test(compute_intersection_13 1)
    define_mpi_test(compute_intersection_22 1)
    define_mpi_test(compute_intersection_23 1)
    
    define_mpi_test(prolongation_12 1)
    define_mpi_test(prolongation_13 1)
    define_mpi_test(prolongation_23 1)
    
    define_mpi_test(intersection_22 1)
    
    define_mpi_test(speed_simple 1)
#     define_mpi_test(speed 1)

