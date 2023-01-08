# makefile for use with GNU make (gmake)
#
# Author: Allan P. Engsig-Karup.
#

include common.mk

.SUFFIXES: .f90 .f

# Source files
# subdirectories - complete directory-tree should be included
#
# Perhaps use perl script sfmakedepend in the future for dependency list?
SOURCES := 
include src/variabledefs/makefile.inc
include src/utilities/makefile.inc
include src/analyticalsolutions/makefile.inc
include src/curvilinear/makefile.inc
include src/fft/makefile.inc
include src/functions/makefile.inc
include src/multigrid/makefile.inc
include src/pressure/makefile.inc
include src/timeintegration/makefile.inc
include src/wrappers/makefile.inc
include src/initialization/makefile.inc
include src/IO/makefile.inc
include src/main/makefile.inc
include src/iterative/makefile.inc
include src/OpenFoam/makefile.inc

# Search paths for source dependencies
VPATH  = src/variabledefs 
VPATH += src/analyticalsolutions
VPATH += src/variabledefs
VPATH += src/analyticalsolutions
VPATH += src/curvilinear
VPATH += src/fft
VPATH += src/functions
VPATH += src/initialization
VPATH += src/IO
VPATH += src/main
VPATH += src/multigrid
VPATH += src/pressure
VPATH += src/timeintegration
VPATH += src/utilities
VPATH += src/wrappers
VPATH += src/iterative
VPATH += src/OpenFoam
VPATH += src/OpenFoam/IO
VPATH += src/OpenFoam/OpenFoam
VPATH += $(BUILDDIR)

# Include directories
INCLUDEDIRS += -I$(BUILDDIR)

# Object files
OBJECTS := $(SOURCES:.f=.o)
OBJECTS := $(OBJECTS:.f90=.o)
#OBJECTS := $(OBJECTS:.f=.o)
OBJECTSNODIR = $(notdir $(OBJECTS))
OBJECTSBUILDDIR = $(addprefix $(BUILDDIR)/,$(OBJECTSNODIR))

# default target creates directory
default:
	@echo "To install OceanWave3D type 'make Debug|Release'"

# Targets for linking
.PHONY: Release
Release: FFLAGS = $(OPTFLAGS)
Release: MODE = "Release"
Release: $(INSTALLDIR)/$(PROGNAME)

.PHONY: Debug
Debug: FFLAGS = $(DBFLAGS)
Debug: MODE = "Debug"
Debug: $(INSTALLDIR)/$(PROGNAME)

$(INSTALLDIR)/$(PROGNAME): $(OBJECTSBUILDDIR) $(INSTALLDIR)
	@if ls *.mod &> /dev/null; then \
	mv -v ./*.mod $(BUILDDIR); \
	cp -v ./thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod $(BUILDDIR); \
	fi
	@echo "*** Starting linking of files for OceanWave3D ($(MODE))... ***"
	@$(FC) $(FFLAGS) -o $(INSTALLDIR)/$(PROGNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS) 	
	@echo "OceanWave3D ($(MODE)) has been built successfully."

.PHONY: all
all: Release

.PHONY: sharedLlib
sharedLib: FFLAGS = $(SHLIBFLAGS)
sharedLib: $(OBJECTSBUILDDIR)
	@if ls *.mod &> /dev/null; then \
	mv -v ./*.mod $(BUILDDIR); \
	cp -v ./thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod $(BUILDDIR); \
	fi
	@echo "*** Starting linking of files for OceanWave3D (Release)... ***"
	@$(FC) $(FFLAGS) -o $(FOAM_USER_LIBBIN)/$(LIBNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS) 	
	@echo "Shared library for OceanWave3D has been built successfully."


# Compile only - compile all source file to build directory
.PHONY: compile
compile: $(OBJECTSBUILDDIR)

.PHONY: clean
clean:
	-rm -r $(BUILDDIR)
	-rm $(INSTALLDIR)/$(PROGNAME)

#
# Special source dependencies
#

# For defining generic rules:
#
# $@ - The current target's full name.
# $* - The current target's file name without a suffix.
# $? - A list of the current target's changed dependencies.
# $< - A single changed dependency of the current target.

# Generic compilation rules
$(BUILDDIR)/%.o: %.f $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR)/%.o: %.f90 $(BUILDDIR)
	@if ls *.mod &> /dev/null; then \
	mv -v *.mod $(BUILDDIR); \
	cp -v thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod $(BUILDDIR); \
	fi
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR)/%.mod: %.f $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR)/%.mod: %.f90 $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR):
	@mkdir -p $@

$(INSTALLDIR):
	@mkdir -p $@
