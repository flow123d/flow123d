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

build/CMakeCache.txt:
	if [ ! -d build ]; then mkdir build; fi
	cd build; cmake ..

cmake: build/CMakeCache.txt

build: cmake
	make -j 4 -C build all

FLOW_BIN=build/bin/flow123d
MPIEXEC_BIN=build/bin/mpiexec

install: build
	if [ -e  $(FLOW_BIN) ]; then rm -f bin/flow123d; cp $(FLOW_BIN) bin; fi
	if [ -e  $(MPIEXEC_BIN) ]; then rm -f bin/mpiexec; cp $(MPIEXEC_BIN) bin; chmod u+x bin/mpiexec; fi

all:  install

# timing of parallel builds (on Core 2 Duo, 4 GB ram)
# N JOBS	O3	g,O0	
# 1 		51s	46s
# 2 		30s	26s
# 4 		31s	27s
# 8 		30s
	
# Remove all generated files
clean: cmake
	make -C build clean

# try to remove all
clean-all: 
	rm -rf build
	make -C third_party clean

# remove everything that is not under version control 
# BE EXTREMELY CAREFUL using this
clean_all_svn:
	bin/svnclean.sh

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

	