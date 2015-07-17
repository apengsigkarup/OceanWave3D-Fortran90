# 
F77 = f77
#F77 = cf77
FFLAGS = -g -Wall 

FILES = rinfo1.o dinfo13.o
## needs library libskit.a in whatever machine version -- 

LIB   = -L/project/darpa/lib/PC -lskit
##LIB   = -L/project/darpa/lib/solaris -lskit

info1.ex: $(FILES) 
	$(F77) -o info1.ex $(FILES) $(LIB) 

clean:
	rm -f *.o *.ex core *.trace *~ 

.f.o  : ;       $(F77) $(FFLAGS) -c $*.f -o $*.o