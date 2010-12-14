#
# $Id$
# $Revision$
# $LastChangedBy$
# $LastChangedDate$
#
# This makefile just provide main rules for: build, documentation and testing
# Build itself takes place in ./src.

include makefile.in

all: bin/mpiexec revnumber
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

revnumber:
	if which "svnversion" ;\
	then echo "#define REVISION \"`svnversion`\"" >include/rev_num.h;\
	else echo "#define REVISION \"`bin/svnversion.sh`SH\"" >include/rev_num.h;\
	fi


clean:
	make -C third_party clean
	make -C src clean
	make -C doc/doxy clean
	rm -f bin/mpiexec

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
