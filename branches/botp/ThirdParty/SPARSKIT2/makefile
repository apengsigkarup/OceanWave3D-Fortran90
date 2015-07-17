#-----------------------------------------------------------------------
#                     S P A R S K I T  makefile 
#-----------------------------------------------------------------------
#
# There are three types of source files in SPARSKIT:
# 1. library source files
# 2. non-library and unsupported source files
# 3. test programs
#
# Simply type "make" to make the SPARSKIT library.
#
# "make all" will additionally make all the unsupported objects.
#
# To make the test programs, use the makefiles in each subdirectory,
# or see the dotests script.
# There are no references to test programs in this makefile.
# Some test programs use the non-library objects from other directories,
# so do a "make all" before making test programs.
#
# Last Updated:  May 9, 1994.

SHELL = /bin/sh
F77 = gfortran
# F77 = f77 
# AR = ar rcv 
AR = ar -rcv 
# OPT = -c 
# OPT = -c -static -64 -mips4 -O3 -r10000
#OPT = -c -O3 -shared -fPIC
OPT = -c -g -shared -fPIC

# Cray options
#F77 = cf77
#AR = bld rv
#OPT = -c -Wf"-dp" 

#
# library objects
#
OBJ =	BLASSM/blassm.o		\
	BLASSM/matvec.o		\
	FORMATS/formats.o	\
	FORMATS/unary.o		\
	INFO/infofun.o		\
	INOUT/inout.o		\
	ITSOL/ilut.o		\
	ITSOL/iters.o		\
	MATGEN/FDIF/genmat.o	\
	MATGEN/FEM/elmtlib2.o	\
	MATGEN/FEM/femgen.o	\
	MATGEN/FEM/meshes.o	\
	MATGEN/MISC/sobel.o	\
	MATGEN/MISC/zlatev.o	\
	ORDERINGS/ccn.o		\
	ORDERINGS/color.o	\
	ORDERINGS/dsepart.o

#
# non-library and unsupported objects
#
OBJ2 =	ITSOL/itaux.o		\
	MATGEN/FDIF/functns.o	\
	MATGEN/FEM/functns2.o	\
	UNSUPP/BLAS1/blas1.o	\
	UNSUPP/MATEXP/exppro.o	\
	UNSUPP/MATEXP/phipro.o	\
	UNSUPP/PLOTS/psgrd.o	\
	UNSUPP/PLOTS/texgrid1.o	\
	UNSUPP/PLOTS/texplt1.o

DIRS =	.			\
	BLASSM			\
	FORMATS			\
	INFO			\
	INOUT			\
	ITSOL			\
	MATGEN/FDIF		\
	MATGEN/FEM		\
	MATGEN/MISC		\
	ORDERINGS		\
	UNSUPP/BLAS1		\
	UNSUPP/MATEXP		\
	UNSUPP/PLOTS

lib libskit.a: $(OBJ) 
	 $(AR) libskit_gfortran.a  $(OBJ) 
##	 ranlib libskit.a

# do not ranlib on some architectures

clean: 
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir; rm -f *~ *.o *.ex .,* fort.* ftn?? *.mat *.pic *.ps *.trace core) ;\
          done

tarit: 
	(cd ..; tar cvf - SPARSKIT2 | gzip -c > SPARSKIT2.tar.gz) 

all: $(OBJ) $(OBJ2) libskit_gfortran_dbg.a

BLASSM/blassm.o: BLASSM/blassm.f
	(cd BLASSM ; $(F77)  $(OPT) blassm.f)
BLASSM/matvec.o: BLASSM/matvec.f
	(cd BLASSM ; $(F77)  $(OPT) matvec.f)
FORMATS/formats.o: FORMATS/formats.f
	(cd FORMATS ; $(F77)  $(OPT) formats.f)
FORMATS/unary.o: FORMATS/unary.f
	(cd FORMATS ; $(F77)  $(OPT) unary.f)
INFO/infofun.o: INFO/infofun.f
	(cd INFO ; $(F77)  $(OPT) infofun.f)
INOUT/inout.o: INOUT/inout.f
	(cd INOUT; $(F77)  $(OPT) inout.f)
ITSOL/ilut.o: ITSOL/ilut.f
	(cd ITSOL; $(F77)  $(OPT) ilut.f)
ITSOL/iters.o: ITSOL/iters.f
	(cd ITSOL; $(F77)  $(OPT) iters.f)
ITSOL/itaux.o: ITSOL/itaux.f
	(cd ITSOL; $(F77)  $(OPT) itaux.f)
MATGEN/FDIF/genmat.o: MATGEN/FDIF/genmat.f
	(cd MATGEN/FDIF ; $(F77)  $(OPT) genmat.f)
MATGEN/FDIF/functns.o: MATGEN/FDIF/functns.f
	(cd MATGEN/FDIF ; $(F77)  $(OPT) functns.f)
MATGEN/FEM/elmtlib2.o: MATGEN/FEM/elmtlib2.f
	(cd MATGEN/FEM ; $(F77)  $(OPT) elmtlib2.f)
MATGEN/FEM/femgen.o: MATGEN/FEM/femgen.f
	(cd MATGEN/FEM ; $(F77)  $(OPT) femgen.f)
MATGEN/FEM/functns2.o : MATGEN/FEM/functns2.f 
	(cd MATGEN/FEM ; $(F77)  $(OPT) functns2.f)
MATGEN/FEM/meshes.o: MATGEN/FEM/meshes.f
	(cd MATGEN/FEM ; $(F77)  $(OPT) meshes.f)
MATGEN/MISC/sobel.o: MATGEN/MISC/sobel.f
	(cd MATGEN/MISC ; $(F77)  $(OPT) sobel.f)
MATGEN/MISC/zlatev.o: MATGEN/MISC/zlatev.f
	(cd MATGEN/MISC ; $(F77)  $(OPT) zlatev.f)
ORDERINGS/ccn.o: ORDERINGS/ccn.f
	(cd ORDERINGS ; $(F77)  $(OPT) ccn.f)
ORDERINGS/color.o: ORDERINGS/color.f
	(cd ORDERINGS ; $(F77)  $(OPT) color.f)
ORDERINGS/dsepart.o: ORDERINGS/dsepart.f
	(cd ORDERINGS ; $(F77)  $(OPT) dsepart.f)
UNSUPP/BLAS1/blas1.o: UNSUPP/BLAS1/blas1.f
	(cd UNSUPP/BLAS1 ; $(F77)  $(OPT) blas1.f)
UNSUPP/MATEXP/exppro.o: UNSUPP/MATEXP/exppro.f
	(cd UNSUPP/MATEXP ; $(F77)  $(OPT) exppro.f)
UNSUPP/MATEXP/phipro.o: UNSUPP/MATEXP/phipro.f
	(cd UNSUPP/MATEXP ; $(F77)  $(OPT) phipro.f)
UNSUPP/PLOTS/psgrd.o : UNSUPP/PLOTS/psgrd.f 
	(cd UNSUPP/PLOTS ; $(F77) $(OPT) psgrd.f)
UNSUPP/PLOTS/texgrid1.o : UNSUPP/PLOTS/texgrid1.f 
	(cd UNSUPP/PLOTS ; $(F77) $(OPT) texgrid1.f)
UNSUPP/PLOTS/texplt1.o : UNSUPP/PLOTS/texplt1.f 
	(cd UNSUPP/PLOTS ; $(F77) $(OPT) texplt1.f)
