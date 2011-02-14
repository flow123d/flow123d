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
# This makefile just provide main rules for: build, documentation and testing
# Build itself takes place in ./src.
#

include makefile.in

all: bin/mpiexec revnumber bin/generic_flow
#	make -C src clean
	make -C third_party all
	make -C src all

	
bin/mpiexec: makefile.in
	# TODO:
	# some time PETSC don't set MPIEXEC well
	# then we should detect existence of PETSC_ARCH/bin/mpiexec
	# last chance is to use system wide mpiexec
	# or die
	 
	if which ${MPIEXEC}; then \
	    echo '#!/bin/bash' > bin/mpiexec; \
	    echo '${MPIEXEC} $$@' >> bin/mpiexec; \
	elif [ -x ${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec ]; then \
	    echo '#!/bin/bash' > bin/mpiexec; \
	    echo '${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec $$@' >> bin/mpiexec; \
	else \
	    echo "Can not guess mpiexec of PETSC configuration"; \
	fi        
	chmod u+x bin/mpiexec

#${BUILD_DIR} default value is "", must be set in makefile.in when running by bitten
bin/generic_flow:
	if [ -z ${MACHINE} ]; then \
		echo "Using default generic_flow"; \
		echo '#!/bin/bash' > bin/generic_flow; \
		echo '${PWD}/${BUILD_DIR}/bin/current_flow.qsub' >> bin/generic_flow; \
	else \
		if [ -e bin/${MACHINE}_flow.sh ]; then \
			echo '#!/bin/bash' > bin/generic_flow; \
			echo '${PWD}/bin/${MACHINE}_flow.sh' >> bin/generic_flow; \
		else \
			echo "make_pbs script for given MACHINE not found, using default"; \
			echo '#!/bin/bash' > bin/generic_flow; \
			echo '${PWD}/bin/current_flow.qsub' >> bin/generic_flow; \
		fi \
	fi
	chmod u+x bin/generic_flow
	
revnumber:
	if which "svnversion" ;\
	then echo "#define REVISION \"`svnversion`\"" >include/rev_num.h;\
	else echo "#define REVISION \"`bin/svnversion.sh`SH\"" >include/rev_num.h;\
	fi

# make package Windows:
# 1) build flow, bcd, ngh, mpiexec
# 2) copy doc, only PDF, wihtout .svn subdirs
# 3) copty tests, only input and output files, no scripts


clean:
	make -C third_party clean
	make -C src clean
	make -C doc/doxy clean
	rm -f bin/mpiexec
	rm -f bin/generic_flow

test: all 
	make -C tests testbase
	
testall: all
	make -C tests testall

online-doc:
	make -C doc/doxy doc
	
%.tst :
	@cd tests;\
	BASE=$*;\
	dir="$${BASE}*";\
	if [ -d $${dir} ];\
	then make -C $${dir} test;\
	else echo "missing test directory $${dir}";\
	fi 
