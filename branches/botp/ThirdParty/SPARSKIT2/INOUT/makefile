FFLAGS = 
F77 = f77

#F77 = cf77
#FFLAGS = -Wf"-dp"

FILES1 = chkio.o
FILES2 = hb2ps.o 
FILES3 = hb2pic.o 

chk.ex: $(FILES1) ../MATGEN/FDIF/functns.o ../libskit.a
	$(F77) $(FFLAGS) -o chk.ex $(FILES1) ../MATGEN/FDIF/functns.o ../libskit.a

hb2ps.ex: $(FILES2) ../libskit32.a
	$(F77) $(FFLAGS) -o hb2ps.ex $(FILES2) ../libskit32.a

hb2pic.ex: $(FILES3) ../libskit.a
	$(F77) $(FFLAGS) -o hb2pic.ex $(FILES3) ../libskit.a

clean:
	rm -f *.o *.ex core *.trace *.pic *.mat *.ps

../MATGEN/FDIF/functns.o:
	(cd ../MATGEN/FDIF; $(F77) $(FFLAGS) -c functns.f)

../libskit.a:
	(cd ..; $(MAKE) $(MAKEFLAGS) libskit.a)

