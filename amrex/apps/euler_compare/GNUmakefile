PRECISION = DOUBLE
PROFILE   = FALSE
DEBUG     = FALSE

DIM       = 2
COMP      = gnu

USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE

AMREX_HOME ?= ../../../external/AMReX
EBASE      := amrex_euler_compare
BL_NO_FORT := TRUE

Bpack := ./Make.package
Blocs := . src include

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

Pdirs := Base Boundary LinearSolvers/MLMG
Ppack += $(foreach dir,$(Pdirs),$(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Bpack)
include $(Ppack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += . src

all: $(executable)
	@echo SUCCESS

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
