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
# $Id: makefile 1102 2011-05-04 14:24:17Z jan.brezina $
# $Revision: 1102 $
# $LastChangedBy: jan.brezina $
# $LastChangedDate: 2011-05-04 16:24:17 +0200 (Wed, 04 May 2011) $
#

include ../../makefile.in
include ../../makefile.include

all: fparser.o foptimizer.o

clean:
	rm -f fparser.o fpoptimizer.o

fparser.o:
	${CPP} ${CFLAGS} -I. -c fparser.cc

foptimizer.o:
	${CPP} ${CFLAGS} -I. -c fpoptimizer.cc