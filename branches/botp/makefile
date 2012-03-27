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

#
#test: $(OBJECTSBUILDDIR)
#	@echo $(OBJECTSBUILDDIR)

# default target creates directory
default:
	-@mkdir -p -v $(BUILDDIR)
	@echo "To install OceanWave3D type 'make Debug|Release'"

# Targets for linking
Release: FFLAGS = $(OPTFLAGS)
Release: $(OBJECTSBUILDDIR)
	@if ls *.mod &> /dev/null; then \
	mv -v ./*.mod $(BUILDDIR); \
	cp -v ./thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod $(BUILDDIR); \
	fi
	@echo "*** Starting linking of files for OceanWave3D (Release)... ***"
	@$(FC) $(FFLAGS) -o $(INSTALLDIR)/$(PROGNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS) 	
	ar -cr ~/lib/libOceanWave3DBuild.a $(BUILDDIR)/*.o
	@echo "OceanWave3D has been built successfully."

Debug: FFLAGS = $(DBFLAGS)
Debug: $(OBJECTSBUILDDIR) 
	@if ls *.mod &> /dev/null; then \
	mv -v *.mod $(BUILDDIR); \
	cp -v thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod $(BUILDDIR); \
	fi
	@echo "*** Starting linking of files for OceanWave3D (Debug)... ***"
	$(FC) $(FFLAGS) -o $(INSTALLDIR)/$(PROGNAME) $(OBJECTSBUILDDIR) $(LIBDIRS) $(LINLIB) $(INCLUDEDIRS) 	
	ar -cr ~/lib/libOceanWave3DBuild.a $(BUILDDIR)/*.o
	@echo "OceanWave3D Debug has been built successfully."

all: Release

# Compile only - compile all source file to build directory
compile: $(OBJECTSBUILDDIR)

clean:
	rm -f $(BUILDDIR)/*.o 
	rm -f $(BUILDDIR)/*.mod 
	rm -f $(INSTALLDIR)/$(PROGNAME)

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
$(BUILDDIR)/%.o: %.f 
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR)/%.o: %.f90
	@if ls *.mod &> /dev/null; then \
	mv -v *.mod $(BUILDDIR); \
	cp -v thirdpartylibs/LIB_VTK_IO/static/lib_vtk_io.mod $(BUILDDIR); \
	fi
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR)/%.mod: %.f
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)

$(BUILDDIR)/%.mod: %.f90
	$(FC) $(FFLAGS) -c $< -o $@ $(INCLUDEDIRS)
