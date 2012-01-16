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
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#

include ../makefile.in

all: armadillo fparser


metis:
	make -C metis-4.0

# use CMAKE_PREFIX_PETH to use our copy of boost
#
# unfortunately I didn't succeed with blas and lapack form PETSC ...
# how to setup armadillo's cmake to use them
armadillo:
	cd armadillo-1.1.92; \
	export PETSC_DIR=${PETSC_DIR}; \
	export PETSC_ARCH=${PETSC_ARCH}; \
	export CMAKE_BLAS_LAPACK_DIR="${PETSC_DIR}/${PETSC_ARCH}/lib"; \
	cmake . -DCMAKE_PREFIX_PATH="../boost-1.42";\
	make

fparser:
	make -C fparser-4.4.3

clean:
	rm -f all
	make -C metis-4.0 realclean
	make -C armadillo-1.1.92 clean
	rm armadillo-1.1.92/CMakeCache.txt
	make -C fparser-4.4.3 clean
	