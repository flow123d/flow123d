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

# This file can be included by makefiles from particular tests
# to use general rules: "test" and "clean"

# disallow default rules (prevents making flow123d from flow123d.sh)
.SUFFIXES:

check_build_tree:
	cd ../.. && bin/git_post_checkout_hook

update: check_build_tree ../../bin/flow123d
	../../bin/tests/run_test.sh ${INI_FILES} ${NPROC} ${FLOW_PARAMS} update

test: check_build_tree ../../bin/flow123d
	../../bin/tests/run_test.sh ${INI_FILES} ${NPROC} ${FLOW_PARAMS}

clean:
	rm -rf output Results test_results; rm -f profiler_*; rm -f vystup.txt
