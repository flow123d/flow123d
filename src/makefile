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
include ../makefile.include

LOCAL_INCLUDE = -I../include -I.

# Important directories
BIN = ../bin
OBJS = ../build/obj
DEPS = ../build/deps

# List of binaries to be compiled
TARGETS = $(BIN)/flow123d

# List of object files of Flow123d
FLOW123D_OBJS = \
$(OBJS)/flow/local_matrix.o\
$(OBJS)/flow/darcy_flow_mh.o\
$(OBJS)/flow/darcy_flow_mh_output.o\
\
$(OBJS)/io/read_ini.o\
$(OBJS)/io/output.o\
$(OBJS)/io/output_vtk.o\
$(OBJS)/io/output_msh.o\
$(OBJS)/io/output_tmp.o\
\
$(OBJS)/mesh/msh_gmshreader.o\
$(OBJS)/mesh/mesh.o\
$(OBJS)/mesh/nodes.o\
$(OBJS)/mesh/elements.o\
$(OBJS)/mesh/sides.o\
$(OBJS)/mesh/edges.o\
$(OBJS)/mesh/boundaries.o\
$(OBJS)/mesh/neighbours.o\
${OBJS}/mesh/topology.o\
\
$(OBJS)/quadrature/quadrature_lib.o\
\
$(OBJS)/reaction/linear_reaction.o\
\
$(OBJS)/semchem/che_semchem.o\
$(OBJS)/semchem/semchem_interface.o\
\
$(OBJS)/system/system.o\
$(OBJS)/system/par_distribution.o\
$(OBJS)/system/sys_function_stack.o\
$(OBJS)/system/sys_profiler.o\
$(OBJS)/system/math_fce.o\
\
$(OBJS)/transport/transport.o\
$(OBJS)/transport/transport_operator_splitting.o\
\
$(OBJS)/equation.o\
$(OBJS)/time_marks.o\
$(OBJS)/materials.o\
$(OBJS)/solve.o\
$(OBJS)/main.o\
$(OBJS)/la_schur.o\
$(OBJS)/xio.o\
$(OBJS)/sparse_graph.o\
$(OBJS)/la_linsys.o\
$(OBJS)/time_governor.o\
$(OBJS)/hc_explicit_sequential.o

#$(OBJS)/convert.o\
#$(OBJS)/mesh/ini_constants_mesh.o\
#$(OBJS)/constantdb.o\
#$(OBJS)/hashes.o\
#$(OBJS)/transport_bcd.o\
#$(OBJS)/concentrations.o\
#$(OBJS)/ppfcs.o\
#$(OBJS)/reaction.o\
#$(OBJS)/reaction.o\
#$(OBJS)/btc.o\

# src subdirs
SUB_DIRS = \
system\
mesh\
io\
quadrature\
semchem\
flow\
transport\
reaction

BUILD_DIRS = $(addprefix $(DEPS)/, $(SUB_DIRS)) $(addprefix $(OBJS)/, $(SUB_DIRS)) 

FLOW123D_DEPS = $(addprefix $(DEPS)/, $(notdir $(FLOW123D_OBJS:.o=.d)))

.PHONY: deps

all: $(TARGETS)

${BUILD_DIRS}:
	echo "make dir:" $@
	mkdir -pf $@ 
	
deps: $(FLOW123D_DEPS)

# General rule to create dependencies files for .c source files
${DEPS}/%.d: %.c ${DEPS}
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(DEPS)/, $(notdir $(@:%.o=%.c)))" $<

# General rule to create dependencies files for .cc source files
${DEPS}/%.d: %.cc ${DEPS}
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(DEPS)/, $(notdir $(@:%.o=%.cc)))" $<

# General rule to create dependencies files for .ccp source files
${DEPS}/%.d: %.ccp  ${DEPS}
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(DEPS)/, $(notdir $(@:%.o=%.ccp)))" $<

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
	rm -f ${TARGETS}
	rm -f ${OBJS}/*.o ${DEPS}/*.d
	rm -f ${OBJS}/flow/*.o
	rm -f ${OBJS}/io/*.o
	rm -f ${OBJS}/mesh/*.o
	rm -f ${OBJS}/reaction/*.o 
	rm -f ${OBJS}/semchem/*.o
	rm -f ${OBJS}/system/*.o
	rm -f ${OBJS}/transport/*.o 
	
