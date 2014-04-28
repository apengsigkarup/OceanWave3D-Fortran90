# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D
LIBNAME  = libOceanWave3D.so

# Installation directory
INSTALLDIR = $(HOME)/bin
LIBINSTALLDIR = $(HOME)/lib


# Build directory where object files are stored 
BUILDDIR = $(PWD)/../build

# The build environment is set either by the choice of a compiler 
# flag, or by creating a block for a specific $USER.  
# Choose the Fortran compiler on this system
# E.g. pathf90, f90, gfortran, gf90, ifort
#FC = gfortran
#FC = gfortran44
#FC = gfortran-4.4
#FC = gf89

USER = apek

# Compiler-dependent section

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

ifeq ($(FC),gf90)
  # hbb home machine
  LIBDIRS  = -L $(HOME)/lib/ 
  LINLIB   = -lharwell -lskit -llapack -lblas
  DBFLAGS  = -pg -g -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none  -fno-automatic
  OPTFLAGS = -O3 -ffree-line-length-none -fno-automatic -ffpe-trap=invalid,zero,overflow
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
  LINLIB   = -lharwell_gfortran -lskit_gfortran -llapack_gfortran -lblas_gfortran 
  DBFLAGS  = -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
endif

# Then this the blocks for specific users (this clobbers the above info.)

ifeq ($(USER),botp)
  # botp kubuntu, 10.04-64bit
  FC       = gfortran-4.6
  LIBDIRS  = -L$(HOME)/OceanWave3D/branches/botp/ThirdParty/lib/ -L/usr/local/include
  LINLIB   = -lharwell_gfortran-4.6 -lskit_gfortran-4.6 -llapack_gfortran-4.6 -lblas -lfftw3
  DBFLAGS  = -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
  OPTFLAGS = -O3 -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
  SHlibFLAGS  = -shared -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
endif

ifeq ($(USER),botpGbar)
  # botp Gbar
  FC       = gfortran
  LIBDIRS  = -L$(HOME)/OceanWave3D/branches/botp/ThirdParty/lib/ 
  LINLIB   = -lharwell_gfortran -lskit_gfortran -llapack_gfortran -lblas
  DBFLAGS  = -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
  SHlibFLAGS  = -shared -O2 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
endif

ifeq ($(USER),apek)
  # MacOS, apek                                                                                                                                                 
  FC       = gfortran
  LIBDIRS  = -L $(HOME)/lib/
# -L/Users/apek/Documents/Fortran/Harwell/lib -L/usr/local/atlas/lib -L/Users/apek/Documents/Fortran/SPARSKIT2                                                    #  LINLIB   = -lharwell -lskit -latlas -llapack -l_VTK_IO                                                                                                         #  LINLIB   = -lharwell -lskit -latlas -l_VTK_IO -framework veclib                                                                                             
  LINLIB   = -lharwell -lskit -latlas -framework veclib
  DBFLAGS  = -pg -g -O -fcheck=all -ffpe-trap=invalid,zero,overflow
  OPTFLAGS = -O2
endif

ifeq ($(USER),olli)
  # olli linux machine
  FC       = gfortran
  LIBDIRS  = -L$(HOME)/lib/ -Ldep/SPARSKIT2/ -Ldep/Harwell/
  LINLIB   = -lharwell -lskit -llapack -lblas 
  DBFLAGS  = -pg -g -O0 -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
  OPTFLAGS = -O3 -ffpe-trap=invalid,zero,overflow -ffree-line-length-none
endif


