#
# $Id: makefile 596 2010-07-21 09:00:40Z jiri.hnidek $
# $Revision: 596 $
# $LastChangedBy: jiri.hnidek $
# $LastChangedDate: 2010-07-21 11:00:40 +0200 (Wed, 21 Jul 2010) $
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
