#!/usr/bin/make
#----------------------------------------------------------------------------------------------------------------------------------
# make init

# shell
SHELL = /bin/bash
# no verbose
$(VERBOSE).SILENT:
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
# COMPILER = gnu   => GNU gfortran (default)
# COMPILER = intel => Intel Fortran
# COMPILER = pgi   => Portland pgf95
# COMPILER = g95   => free g95
#
# SHARED   = yes/no
# DEBUG    = yes/no
# F03STD   = yes/no
# OPTIMIZE = yes/no
# OPENMP   = yes/no
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
SHARED   = no
COMPILER = gnu
DEBUG    = no
F03STD   = no
OPTIMIZE = no
OPENMP   = no
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
DSRC  = ./source/
DOBJ  = ./obj/
DMOD  = ./mod/
VPATH = $(DSRC) $(DOBJ) $(DMOD)
ifeq "$(SHARED)" "yes"
  DLIB    = ./shared/
  MAIN    = $(DLIB)Lib_VTK_IO.so
  MAKELIB = $(FC) $(OPTSL) $(DOBJ)lib_vtk_io.o -o $(MAIN)
else
  DLIB    = ./static/
  MAIN    = $(DLIB)Lib_VTK_IO.a
  MAKELIB = ar -rcs $(MAIN) $(DOBJ)lib_vtk_io.o ; ranlib $(MAIN)
endif

MKDIRS = $(DOBJ) $(DMOD) $(DLIB)
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
ifeq "$(COMPILER)" "gnu"
  FC    = gfortran
  OPTSC = -cpp -c -Dr16p -J$(DMOD)
  OPTSL =
  ifeq "$(SHARED)" "yes"
    OPTSC := $(OPTSC) -fPIC
    OPTSL := $(OPTSL) -shared
  endif
  ifeq "$(DEBUG)" "yes"
    OPTSC := $(OPTSC) -O0 -Wall -Warray-bounds -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow,underflow,precision,denormal
    OPTSL := $(OPTSL) -O0 -Wall -Warray-bounds -fcheck=all -fbacktrace -ffpe-trap=invalid,overflow,underflow,precision,denormal
  endif
  ifeq "$(F03STD)" "yes"
    OPTSC := $(OPTSC) -std=f2008 -fall-intrinsics
    OPTSL := $(OPTSL) -std=f2008 -fall-intrinsics
  endif
  ifeq "$(OPTIMIZE)" "yes"
    OPTSC := $(OPTSC) -O3
    OPTSL := $(OPTSL) -O3
  endif
  ifeq "$(OPENMP)" "yes"
    OPTSC := $(OPTSC) -fopenmp
    OPTSL := $(OPTSL) -fopenmp
  endif
endif
ifeq "$(COMPILER)" "intel"
  FC    = ifort
  OPTSC = -cpp -c -Dr16p -module $(DMOD)
  OPTSL =
  ifeq "$(SHARED)" "yes"
    OPTSC := $(OPTSC) -fpic
    OPTSL := $(OPTSL) -shared
  endif
  ifeq "$(DEBUG)" "yes"
    CHK = -check all -check noarg_temp_created
    DEB = -debug all
    WRN = -warn all
    OPTSC := $(OPTSC) -O0 -fpe-all=0 -fp-stack-check -traceback $(WRN) $(CHK) $(DEB)
    OPTSL := $(OPTSL) -O0 -fpe-all=0 -fp-stack-check -traceback $(WRN) $(CHK) $(DEB)
  endif
  ifeq "$(F03STD)" "yes"
    OPTSC := $(OPTSC) -std03
    OPTSL := $(OPTSL) -std03
  endif
  ifeq "$(OPTIMIZE)" "yes"
    OPTSC := $(OPTSC) -O3 -ipo
    OPTSL := $(OPTSL) -O3 -ipo
  endif
  ifeq "$(OPENMP)" "yes"
    OPTSC := $(OPTSC) -openmp
    OPTSL := $(OPTSL) -openmp
  endif
endif
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
VTK_IO: $(MKDIRS) $(DOBJ)lib_vtk_io.o
	@$(MAKELIB)
	@mv -f $(DMOD)lib_vtk_io.mod $(DLIB)
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
COTEXT = "Compiling $<"
$(DOBJ)lib_vtk_io.o: Lib_VTK_IO.f90 $(DOBJ)ir_precision.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) $< -o $@

$(DOBJ)ir_precision.o: IR_Precision.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) $< -o $@
#----------------------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------------------------------
clean:
	@echo "deleting objs and mods"
	@rm -rf $(DOBJ) $(DMOD)

cleanall: clean
	@echo "deleting libs"
	@rm -fr $(DLIB) ./shared/
#----------------------------------------------------------------------------------------------------------------------------------
