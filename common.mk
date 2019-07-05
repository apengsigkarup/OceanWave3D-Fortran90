# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
<<<<<<< HEAD
PROGNAME = OceanWave3D-Fortran
LIBNAME  = libOceanWave3D_botp.so
=======
PROGNAME = OceanWave3D
LIBNAME  = libOceanWave3D.so
>>>>>>> bb7d99eaba58f90227126d8f4f705ca8218bbfbb

# Installation directory
INSTALLDIR = $(HOME)/bin
LIBINSTALLDIR = $(HOME)/lib

# Build directory where object files are stored 
<<<<<<< HEAD
BUILDDIR = $(PWD)/buildDevelop
=======
BUILDDIR = $(PWD)/build

# The build environment is set either by the choice of a compiler 
# flag, or by creating a block for a specific $USER.  
# Choose the Fortran compiler on this system
# E.g. pathf90, f90, gfortran, gf90, ifort
#FC = gfortran
#FC = gfortran44
#FC = gfortran-4.4
#FC = gf90

USER = botp-dev

# First the blocks based on compiler name:  

ifeq ($(FC),gfortran44)
  # gbar linux machines
  LIBDIRS  = -L $(HOME)/lib/ 
  LINLIB   = -lharwell -lskit -llapack -lblas
  DBFLAGS  = -pg -g -O0 -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
endif

ifeq ($(FC),ifort)
  # hbb work machine with intel compiler
  LIBDIRS  = -L/usr/local/lib/ 
  LINLIB   = -lharwell_intel -lskit_intel -llapack_intel -lblas_intel
  DBFLAGS  = -g -CB -fpe0 -fpstkchk  -traceback 
  OPTFLAGS = -O3 -tpp7
#  OPTFLAGS = -g -fast
endif


ifeq ($(FC),gfortran)
  # MacOS, apek
  LIBDIRS  = -L $(HOME)/lib/
# -L/Users/apek/Documents/Fortran/Harwell/lib -L/usr/local/atlas/lib -L/Users/apek/Documents/Fortran/SPARSKIT2
#  LINLIB   = -lharwell -lskit -latlas -llapack -l_VTK_IO
#  LINLIB   = -lharwell -lskit -latlas -l_VTK_IO -framework veclib
  LINLIB   = -lharwell -lskit -latlas -framework veclib
  DBFLAGS  = -pg -g -O -fcheck=all -ffpe-trap=invalid,zero,overflow
  OPTFLAGS = -O2
endif

ifeq ($(FC),pathf90)
  # Niflheim cluster, DTU, apek
  LIBDIRS  = -L/opt/acml3.5.0/pathscale64/lib -L$(HOME)/lib
  LINLIB   = -lharwellFPATH90 -lskitPATHF90 -llinpackPATHF90 -lacml
  DBFLAGS  = -pg -g -O -static -woffoptions
  OPTFLAGS = -O2 -static -woffoptions -ffortran-bounds-check
endif

ifeq ($(FC),f90)
  # gbar, DTU, apek
  LIBDIRS  = -L$(HOME)/lib/ 
  LINLIB   = -lharwell -lskit -xlic_lib=sunperf
  DBFLAGS  = -pg -g -O0 # -static -woffoptions
  OPTFLAGS = -O -fast 
#  OPTFLAGS = -g -fast
endif
>>>>>>> bb7d99eaba58f90227126d8f4f705ca8218bbfbb

FC = gfortran
PLATFORM = gbar
# Then the blocks for specific users (this clobbers the above info.)
# -lma41  -lma27 -lmi24 -lhsl_mi20
ifeq ($(PLATFORM),gbar)
  LINLIB = 
  DBFLAGS = -I/usr/lib64/gfortran/modules
  OPTFLAGS = -I/usr/lib64/gfortran/modules
endif
<<<<<<< HEAD
ifeq ($(PLATFORM),jess)
  LINLIB = -lma41
  DBFLAGS = -I$(EBROOTHDF5)/include
  OPTFLAGS = -I$(EBROOTHDF5)/include
=======

ifeq ($(USER),botp-dev)
  # botp kubuntu, 10.04-64bit
  FC       = gfortran
  LIBDIRS  = -L$(HOME)/lib/ 
  LINLIB   = -lskit_gfortran -ltmglib_gfortran -llapack_gfortran -lblas 
  DBFLAGS  = -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
  SHLIBFLAGS  = -shared -O2 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
>>>>>>> bb7d99eaba58f90227126d8f4f705ca8218bbfbb
endif

ifeq ($(FC),gfortran)
  # fabpi machine, gfortran
  LIBDIRS  = -L$(HOME)/lib/  -L$(MODULE_OPENBLAS_LIB_DIR)
  LINLIB   += -lharwell -lskit -lopenblas -lhdf5 -lhdf5_fortran -lhdf5_hl -lma41
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
