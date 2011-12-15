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
include makefile.include


all: bin/mpiexec revnumber bin/current_flow
	make -C third_party all
	make -j 4 -C src all
	
bin/mpiexec: makefile.in
	# TODO:
	# some time PETSC don't set MPIEXEC well
	# then we should detect existence of PETSC_ARCH/bin/mpiexec
	# last chance is to use system wide mpiexec
	# or die
	 
	if which "${MPIEXEC}"; then \
	    echo '#!/bin/bash' > bin/mpiexec; \
	    echo '"${MPIEXEC}" "$$@"' >> bin/mpiexec; \
	elif [ -x "${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec" ]; then \
	    echo '#!/bin/bash' > bin/mpiexec; \
	    echo '"${PETSC_DIR}/${PETSC_ARCH}/bin/mpiexec" "$$@"' >> bin/mpiexec; \
	else \
	    echo "Can not guess mpiexec of PETSC configuration"; \
	fi        
	chmod u+x bin/mpiexec

bin/current_flow:
	if [ -z "${MACHINE}" ]; then \
		echo "Using default: current_flow"; \
		echo '#!/bin/bash' > bin/current_flow; \
		echo "\"`pwd`/bin/generic_flow.sh\"" >> bin/current_flow; \
	else \
		if [ -e "bin/stub/${MACHINE}_flow.sh" ]; then \
			echo '#!/bin/bash' > bin/current_flow; \
			echo '"`pwd`/bin/stub/${MACHINE}_flow.sh"' >> bin/current_flow; \
		else \
			echo "script for given MACHINE not found, using default"; \
			echo '#!/bin/bash' > bin/current_flow; \
			echo '"`pwd`/bin/stub/generic_flow.sh"' >> bin/current_flow; \
		fi \
	fi
	chmod u+x bin/current_flow
		#echo '"${PWD}/${BUILD_DIR}/bin/generic_flow.sh"' >> bin/current_flow; \
	
revnumber:
	if which "svnversion" ;\
	then echo "#define REVISION \"`svnversion`\"" >__tmp__ren_num.h;\
	else echo "#define REVISION \"`bin/svnversion.sh`SH\"" >__tmp__ren_num.h;\
	fi ;\
	if test -r "include/rev_num.h" && diff "__tmp__ren_num.h" "include/rev_num.h" >/dev/null; \
	then rm __tmp__ren_num.h; \
	else mv -f __tmp__ren_num.h include/rev_num.h; \
	fi
	


# Remove all generated files
clean:
	make -C src clean
	make -C doc/doxy clean
	make -C tests clean
	rm -f bin/mpiexec
	rm -f bin/current_flow

# Make all tests	
testall:
	make -C tests testall

# Make only certain test (eg: make 01.tst will make first test)
%.tst :
	make -C tests $*.tst

# Create doxygen documentation
online-doc:
	make -C doc/doxy doc

clean_tests:
	make -C tests clean

ngh:
	make -C bin/ngh all

bcd:
	make -C bin/bcd all

clean_util:
	make -C bin/bcd clean
	make -C bin/ngh clean

lbuild=linux_build
linux_package: clean clean_tests clean_util all bcd ngh
	# copy bin
	rm -rf $(lbuild)
	mkdir -p $(lbuild)/bin/mpich
	mpiexec=`cat bin/mpiexec |grep mpiexec |sed 's/ ".*$$//'|sed 's/"//g'`;\
	cp "$${mpiexec}" $(lbuild)/bin/mpich/mpiexec
	cp -r bin/flow123d bin/flow123d.sh bin/ndiff bin/tests bin/ngh/bin/ngh bin/bcd/bin/bcd $(lbuild)/bin
	cp -r bin/paraview $(lbuild)/bin
	# copy doc
	mkdir $(lbuild)/doc
	cp -r doc/articles doc/reference_manual/flow123d_doc.pdf doc/petsc_options_help $(lbuild)/doc
	mkdir $(lbuild)/doc/ngh
	mkdir $(lbuild)/doc/bcd
	cp bin/ngh/doc/* $(lbuild)/doc/ngh
	cp bin/bcd/doc/* $(lbuild)/doc/bcd
	# copy tests
	cp -r tests $(lbuild)

linux_pack:
	cd $(lbuild); tar -cvzf ../flow_build.tar.gz .

	