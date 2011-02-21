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

# Include global settings
include ../makefile.in

LOCAL_INCLUDE = -I../include

# Important directories
BIN = ../bin
OBJS = ../build/obj
DEPS = ../build/deps

# List of binaries to be compiled
TARGETS = $(BIN)/flow123d

# List of object files of Flow123d
FLOW123D_OBJS = \
$(OBJS)/mesh/ini_constants_mesh.o\
$(OBJS)/transport_bcd.o\
$(OBJS)/reaction.o\
$(OBJS)/system.o\
$(OBJS)/read_ini.o\
$(OBJS)/constantdb.o\
$(OBJS)/math_fce.o\
$(OBJS)/msh_gmshreader.o\
$(OBJS)/nodes.o\
$(OBJS)/elements.o\
$(OBJS)/sides.o\
$(OBJS)/edges.o\
$(OBJS)/hashes.o\
$(OBJS)/boundaries.o\
$(OBJS)/neighbours.o\
$(OBJS)/mesh.o\
$(OBJS)/materials.o\
$(OBJS)/local_matrix.o\
$(OBJS)/darcy_flow_mh.o\
$(OBJS)/output.o\
$(OBJS)/concentrations.o\
$(OBJS)/darcy_flow_mh_output.o\
$(OBJS)/solve.o\
$(OBJS)/problem.o\
$(OBJS)/main.o\
$(OBJS)/btc.o\
$(OBJS)/par_distribution.o\
$(OBJS)/la_schur.o\
$(OBJS)/xio.o\
$(OBJS)/sparse_graph.o\
$(OBJS)/convert.o\
$(OBJS)/sys_function_stack.o\
$(OBJS)/transport.o\
$(OBJS)/reaction.o\
$(OBJS)/profiler.o\
$(OBJS)/la_linsys.o\
$(OBJS)/time_governor.o\
$(OBJS)/semchem/che_semchem.o\
$(OBJS)/semchem/interfaceN.o\
$(OBJS)/semchem/Cti_ichNew.o

# $(OBJS)/ppfcs.o\

FLOW123D_DEPS = $(addprefix $(DEPS)/, $(notdir $(FLOW123D_OBJS:.o=.d)))

.PHONY: deps


all: $(TARGETS)

deps: $(FLOW123D_DEPS)

# General rule to create dependencies files for .c source files
${DEPS}/%.d: %.c
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# General rule to create dependencies files for .cc source files
${DEPS}/%.d: %.cc
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# General rule to create dependencies files for .ccp source files
${DEPS}/%.d: %.ccp
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# Include generated files with dependencies
-include $(FLOW123D_DEPS)

# General rule to compile .c source files (pure C)
${OBJS}/%.o: %.c ../makefile.in
	${CPP} -c -o $@ ${CFLAGS} ${CDEFS} ${LOCAL_INCLUDE} ${INCLUDE} $<

# General rule to compile .cc source files (mixture of C and C++)
${OBJS}/%.o: %.cc ../makefile.in
	${CPP} -c -o $@ ${CFLAGS} ${CDEFS} ${LOCAL_INCLUDE} ${INCLUDE} $<	

# General rule to compile .cpp source files (pure C++)
${OBJS}/%.o: %.cpp ../makefile.in
	${CPP} -c -o $@ ${CFLAGS} ${CDEFS} ${LOCAL_INCLUDE} ${INCLUDE} $<	

# Rule to compile one of the target
${BIN}/flow123d: ${FLOW123D_OBJS}
	${CPP_LINKER} ${CFLAGS} -o $@ $^ ${LIBS}

# The rule to delete all generated files	
clean: 
	rm -f ${OBJS}/*.o ${DEPS}/*.d ${TARGETS}
