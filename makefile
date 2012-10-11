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

all:  build_all install

FLOW_BIN=flow123d
INTERPOLATE_BIN=interpolation
MPIEXEC_BIN=mpiexec

ifndef N_JOBS
  N_JOBS=4
endif  

# install all binaries form build tree to './bin' dir
install: 
	if [ -e  "build/$(INTERPOLATE_BIN)" ]; then rm -f bin/$(INTERPOLATE_BIN); cp "build/$(INTERPOLATE_BIN)" bin; fi
	if [ -e  "build/$(FLOW_BIN)" ]; then rm -f bin/$(FLOW_BIN); cp "build/$(FLOW_BIN)" bin; fi
	if [ -e  "build/$(MPIEXEC_BIN)" ]; then rm -f bin/$(MPIEXEC_BIN); cp "build/$(MPIEXEC_BIN)" bin; chmod a+x bin/mpiexec; fi


# run first cmake
build/CMakeCache.txt:
	if [ ! -d build ]; then mkdir build; fi
	cd build; cmake ..

# This target builds links in directory test_units and its subdirectories 
# to generated makefiles in the build directory. 
# This way we can run tests from the source tree and do not have problems with deleted
# current directory in shell if we are forced to use make clean-all.
create_unit_test_links:
	for f in  `find test_units/ -name CMakeLists.txt`; do ln -sf "$${PWD}/build/$${f%/*}/Makefile" "$${f%/*}/makefile";done

# This target only configure the build process.
# Useful for building unit tests without actually build whole program.
cmake: build/CMakeCache.txt  create_unit_test_links



build_all: build_flow123d  ngh

flow123d:  build_flow123d  install


# timing of parallel builds (on Core 2 Duo, 4 GB ram)
# N JOBS	O3	g,O0	
# 1 		51s	46s
# 2 		30s	26s
# 4 		31s	27s
# 8 		30s
build_flow123d: cmake
	make -j $(N_JOBS) -C build flow123d

	
interpolation: build_interpolation install
	
build_interpolation: 
	make -j $(N_JOBS) -C build interpolation


	
# Remove all generated files
clean: cmake
	make -C build clean

# try to remove all
clean-all: 
	rm -f bin/${FLOW_BIN}
	rm -f bin/${MPIEXEC_BIN}
	rm -f bin/${INTERPOLATE_BIN}
	rm -rf build
	for f in  `find test_units/ -name makefile`; do rm -f "$${f}";done
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

	
