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

# This makefile defines basic rules for 
# compiling unit tests in their direcotries. In particular directories are relative to the directory of a unit test


# Include global settings
include ../../makefile.in

LOCAL_INCLUDE = -I../../include

# Important directories

# binery of a unit test
BIN = .
# flow123d objs
OBJS = ../../build/obj
# flow123d dependency files
DEPS = ../../build/deps
# flow123d sources
SRC=../../src

# List of binaries to be compiled
TARGETS = $(BIN)/test

# List of object files needed by a particular unit test
#TEST_OBJS =

TEST_DEPS = $(addprefix $(DEPS)/, $(notdir $(TEST_OBJS:.o=.d)))

.PHONY: deps

deps: $(TEST_DEPS)
	echo " all deps"
	echo ${TEST_DEPS}

# General rule to create dependencies files for .cc source files
${DEPS}/%.d: ${SRC}/%.cc
	echo "generating dependency file " $@
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# General rule to create dependencies files for .ccp source files
${DEPS}/%.d: ${SRC}/%.ccp
	echo "generating dependency file " $@
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# Include generated files with dependencies
-include $(TEST_DEPS)

# General rule to compile .cc source files (mixture of C and C++)
${OBJS}/%.o: ${SRC}/%.cc
	${CPP} -c -o $@ ${CFLAGS} ${CDEFS} ${LOCAL_INCLUDE} ${INCLUDE} $<	

# General rule to compile .cpp source files (pure C++)
${OBJS}/%.o: ${SRC}/%.cpp
	${CPP} -c -o $@ ${CFLAGS} ${CDEFS} ${LOCAL_INCLUDE} ${INCLUDE} $<	

# Rule to compile one of the target
test: ${TEST_OBJS} 
	${CPP} -c -o test.o ${DEFS} ${LOCAL_INCLUDE} ${INCLUDE} test.cpp
	${CPP_LINKER} ${CFLAGS} -o test test.o ${TEST_OBJS} ${LIBS}

# The rule to delete all generated files	
clean: 
	rm -f ${OBJS}/*.o ${DEPS}/*.d ${TARGETS} test.o
