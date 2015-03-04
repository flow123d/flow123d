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

# Make all tests
test-all: 02.tst 03.tst 05.tst 06.tst 08.tst 10.tst 11.tst 12.tst 13.tst 14.tst 15.tst 16.tst 17.tst 18.tst 19.tst 20.tst 21.tst 22.tst 23.tst 24.tst

# Rule to clean generated files in test directories
clean:
	for d in * ;\
	do \
	if [ -d "$${d}" ];\
	then \
		make -C "$${d}" clean;\
	fi;\
	done
	

# General rule to run test given only by number
%.tst:
	BASE=$*;DIR="`pwd`/`echo $${BASE}_*`";if test -d "$${DIR}";then make -C "$${DIR}" test;else exit 1;fi
