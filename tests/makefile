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

# this file can be included by makefiles from particular tests
# to use general test: and clean: rules

### rules for tests
test:
	../../bin/tests/run_test.sh ${INI_FILES} ${NPROC} ${FLOW_PARAMS}

clean:
clean:
	rm -f profiler_*
	rm -f output/*

### global rules

# Make all tests
testall: 01.tst 02.tst


# Rules for certain tests
#01.tst:
#	make -C 01_steady_flow_123d test
#02.tst:
#	make -C 02_transport_12d test

# general rule to run test given only by number
%.tst:
	BASE=$*;\
	dir="`pwd`/`echo $${BASE}*`";\
	echo $${dir};\
	if test -d "$${dir}";\
	then make -C "$${dir}" test;\
	else echo "missing test directory $${dir}";\
	exit 1;\
	fi


# Rule to clean generated files in test directories
cleanall:
	make -C 01_steady_flow_123d clean
	make -C 02_transport_12d clean
