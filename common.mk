# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D-Fortran-Release
LIBNAME  = libOceanWave3D_botp.so

# Installation directory
INSTALLDIR = $(HOME)/bin
LIBINSTALLDIR = $(HOME)/lib

# Build directory where object files are stored 
BUILDDIR = $(PWD)/buildRelease

FC = gfortran
PLATFORM = jess
# Then the blocks for specific users (this clobbers the above info.)
# -lma41  -lma27 -lmi24 -lhsl_mi20
ifeq ($(PLATFORM),gbar)
  LINLIB = 
  DBFLAGS = -I/usr/lib64/gfortran/modules
  OPTFLAGS = -I/usr/lib64/gfortran/modules
endif
ifeq ($(PLATFORM),jess)
  LINLIB = -lma41
  DBFLAGS = -I$(EBROOTHDF5)/include -gdwarf-2
  OPTFLAGS = -I$(EBROOTHDF5)/include
endif
   
ifeq ($(FC),gfortran)
  # fabpi machine, gfortran
  LIBDIRS  = -L$(HOME)/lib/  -L$(MODULE_OPENBLAS_LIB_DIR)
  # LINLIB   += -lharwell -lskit -lopenblas -lhdf5 -lhdf5_fortran -lhdf5_hl
  LINLIB   += -lskit -lopenblas -lhdf5 -lhdf5_fortran -lhdf5_hl
  DBFLAGS  += -pg -g -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none  -fno-automatic
  OPTFLAGS += -pg -O3 -ffree-line-length-none -fno-automatic -ffpe-trap=invalid,zero,overflow 
endif

ifeq ($(FC),ifort)
  # hbb work machine with intel compiler
  LIBDIRS  = -L$(HOME)/lib/  -L$(MODULE_OPENBLAS_LIB_DIR)
  #  LINLIB   = -mkl -lharwell_intel -lskit_intel -llapack_intel -lblas_intel
  LINLIB   += -lopenblas -lma27 -lma41 
  DBFLAGS  += -g -CB -fpe0 -fpstkchk  -traceback 
  OPTFLAGS += -O3 -tpp7
endif

print-%  : ; @echo $* = $($*)
