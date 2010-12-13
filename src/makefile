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
$(OBJS)/transport.o\
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
$(OBJS)/initials.o\
$(OBJS)/materials.o\
$(OBJS)/local_matrix.o\
$(OBJS)/calculation.o\
$(OBJS)/darcy_flow_mh.o\
$(OBJS)/output.o\
$(OBJS)/concentrations.o\
$(OBJS)/postprocess.o\
$(OBJS)/preprocess.o\
$(OBJS)/solve.o\
$(OBJS)/problem.o\
$(OBJS)/main.o\
$(OBJS)/ppfcs.o\
$(OBJS)/btc.o\
$(OBJS)/reaction.o\
$(OBJS)/par_distribution.o\
$(OBJS)/la_schur.o\
$(OBJS)/xio.o\
$(OBJS)/sparse_graph.o\
$(OBJS)/convert.o\
$(OBJS)/sys_function_stack.o\
$(OBJS)/la_linsys.o

FLOW123D_DEPS = $(addprefix $(DEPS)/, $(notdir $(FLOW123D_OBJS:.o=.d)))

.PHONY: deps

all: $(TARGETS)

deps: $(FLOW123D_DEPS)

# General rule to create dependencies files for .cc source files
${DEPS}/%.d: %.cc
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# General rule to create dependencies files for .ccp source files
${DEPS}/%.d: %.ccp
	${CPP} -o $@ ${LOCAL_INCLUDE} ${INCLUDE} -MM -MG -MP -MT "$(addprefix $(OBJS)/, $(notdir $(@:%.d=%.o)))" $<

# Include generated files with dependencies
-include $(FLOW123D_DEPS)

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
