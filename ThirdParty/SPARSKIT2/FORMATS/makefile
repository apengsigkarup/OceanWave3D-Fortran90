FFLAGS =
F77 = f77

#F77 = cf77
#FFLAGS = -Wf"-dp"

FILES1 = chkun.o 
FILES2 = chkfmt1.o 
FILES3 = rvbr.o

un.ex:	$(FILES1) ../MATGEN/FDIF/functns.o ../libskit.a 
	$(F77) $(FFLAGS) -o un.ex $(FILES1) ../MATGEN/FDIF/functns.o ../libskit.a 

fmt.ex:	$(FILES2) ../MATGEN/FDIF/functns.o ../libskit.a 
	$(F77) $(FFLAGS) -o fmt.ex $(FILES2) ../MATGEN/FDIF/functns.o ../libskit.a 

rvbr.ex: $(FILES3) ../MATGEN/FDIF/functns.o ../libskit.a 
	$(F77) $(FFLAGS) -o rvbr.ex $(FILES3) ../MATGEN/FDIF/functns.o ../libskit.a

clean:
	rm -f *.o *.ex core *.trace

../MATGEN/FDIF/functns.o:
	(cd ../MATGEN/FDIF; $(F77) $(FFLAGS) -c functns.f)

../libskit.a:
	(cd ..; $(MAKE) $(MAKEFLAGS) libskit.a)

