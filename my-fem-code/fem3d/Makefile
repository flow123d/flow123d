lib = ad.o \
      mesh3d.o \
      integration3d.o \
      solution.o \
      sparse.o \
      stateproblem.o \
      wmanager.o

targets = ../lib/$(TARGET)/libfem3d.a

default: fem

TOPDIR = ../
include ../Make3d.inc.$(PLATFORM)



fem: $(lib)
	rm -f $(targets)
	ar cr $(targets) $(lib)
	ranlib $(targets)

