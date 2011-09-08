targets = stokes3d stokes2d

TOPDIR = ../
include ../Make3d.inc.$(PLATFORM)

BUILDLIBS = $(FEM) $(SLU) $(UMFPACK) $(BLAS) $(MUPARSER) $(VTK) -lncurses

$(TARGET)/stokes2d: $(TARGET)/stokes2d.o $(BUILDLIBS)

$(TARGET)/stokes3d: $(TARGET)/stokes3d.o $(BUILDLIBS)
