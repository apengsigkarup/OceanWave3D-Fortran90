# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D-Fortran
LIBNAME  = libOceanWave3D_botp.so

# Installation directory
INSTALLDIR = $(HOME)/bin
LIBINSTALLDIR = $(HOME)/lib

# Build directory where object files are stored 
BUILDDIR = $(PWD)/build

FC = gfortran

# Then the blocks for specific users (this clobbers the above info.)
ifeq ($(FC),gfortran)
  # fabpi machine, gfortran
  FC=gfortran
  # LIBDIRS  = -L$(HOME)/lib/ 
  LIBDIRS  = -L$(HOME)/lib/  -L$(MODULE_OPENBLAS_LIB_DIR)
  LINLIB   = -lharwell -lskit -lopenblas -lma41  -lma27 -lmi24 -lhsl_mi20 
  DBFLAGS  = -pg -g -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none  -fno-automatic
  OPTFLAGS = -pg -O3 -ffree-line-length-none -fno-automatic -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(FC),ifort)
  # hbb work machine with intel compiler
  LIBDIRS  = -L$(HOME)/lib/  -L$(MODULE_OPENBLAS_LIB_DIR)
#  LINLIB   = -mkl -lharwell_intel -lskit_intel -llapack_intel -lblas_intel
  LINLIB   = -lopenblas -lma27 -lma41 
  DBFLAGS  = -g -CB -fpe0 -fpstkchk  -traceback 
  OPTFLAGS = -O3 -tpp7
endif
