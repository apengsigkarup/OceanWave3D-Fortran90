#-------------
#-- Flags
#-------------

#LINUX
#FF = /usr/absoft/bin/f90
#FFLAGS = -O -B100
#installd = /usr/local/lib

#HPC SUN machines
#FF = f90
#FFLAGS = -C -xO3 -xchip=ultra3 
#installd = $(HOME)/lib

#Kubuntu 12.04, botp
FF = gfortran
FFLAGS = -O3 -ffree-line-length-none -pedantic -fno-automatic -shared -fPIC
installd = ../lib/

#IBM1
#FF = f90
#FFLAGS = -O3  -Q -qarch=pwr3 -qtune=pwr3 -qflttrap=zerodivide:overflow:invalid:enable 
#FFLAGS = -c -g -C -qflttrap=zerodivide:overflow:invalid:enable 

#installd = /usr/local/lib

#-------------
#-- Rules
#-------------
.f.o:
	$(FF) -c $(FFLAGS) $*.f

#-------------
#-- what
#-------------
target = libharwell_$(FF)


#------------
#-- objects
#------------
objects = MA30AD.o MC13D.o MC20AD.o MA28AD.o MC21A.o MC22AD.o MC23AD.o MC24AD.o MA48.o MA41.o MA50.o MC21.o MA13.o FD05.o MC47.o MC41.o MC49.o MC51.o MI11.o MI26.o MA41SP.o MA41depSP.o FD05SP.o ID05SP.o ZA02SP.o MC51A.o MI24.o MI24Test.o

#------------
#-- targets
#------------
$(target):  $(objects) 
	ar -rc $(installd)/$(target).a $(objects)

#rm $(installd)/$(target).a
#ar -rc $(installd)/$(target).a $(objects)

clean: 
	rm $(objects) 

