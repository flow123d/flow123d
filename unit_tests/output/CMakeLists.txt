# 
# Copyright (C) 2007 Technical University of Liberec.  All rights reserved.
#
# Please make a following refer to Flow123d on your project site if you use the program for any purpose,
# especially for academic research:
# Flow123d, Research Center: Advanced Remedial Technologies, Technical University of Liberec, Czech Republic
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


set(libs coupling_lib mesh_lib flow123d_system)
add_test_directory("${libs}")

define_mpi_test( reader_cache 1 )

define_mpi_test( output 1 )
define_mpi_test( output_vtk 1)
define_mpi_test( output_msh 1)
define_mpi_test( output_mesh 1)
define_mpi_test( observe 1)
define_mpi_test( observe 2)
define_mpi_test( parallel_output 2 )
define_mpi_test( element_data_cache 1 )
define_mpi_test( element_data_cache 2 )
define_mpi_test( element_data_cache 3 60)
define_mpi_test( optimizer_output 1 )
