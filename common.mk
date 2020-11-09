# For inclusion in makefile for use with GNU make (gmake)
# 
# Purpose: Modify to local configuration by user.
#

# Program name
PROGNAME = OceanWave3D
LIBNAME  = libOceanWave3D.so

# Installation directory
INSTALLDIR ?= install/bin

# Build directory where object files are stored 
BUILDDIR = build

FC = gfortran
LIBDIRS  ?= -Llib
LINLIB   ?= -llapack  -lskit -lblas -lharwell
DBFLAGS  ?= -pg -g -O0 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none 
OPTFLAGS ?= -O3 -fPIC -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all
SHLIBFLAGS  ?= -shared -O2 -fPIC -fbounds-check -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -fstack-protector-all

