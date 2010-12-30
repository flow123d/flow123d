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

copy:
	cp -r ../../../tests/test1 ./
	cp -r ../../../tests/test2 ./
	cp -r ../../../tests/test3 ./
	cp -r ../../../tests/test4 ./

testbase:
	make -C 01_steady_flow_123d test	

testall: 
	make -C 01_steady_flow_123d test
	make -C 02_transport_12d test
	make -C 03_transport_small_12d test
	make -C 04_density_driven_12d test

clean:
	rm -rf ./test1
	rm -rf ./test2
	rm -rf ./test3
	rm -rf ./test4

	
%.tst:
	@BASE=$*;\
	dir="$${BASE}*";\
	if [ -d $${dir} ];\
	then make -C $${dir} test;\
	else echo "missing test directory $${dir}";\
	fi