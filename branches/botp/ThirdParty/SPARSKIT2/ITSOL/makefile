FFLAGS =
F77 = f77

#F77 = cf77
#FFLAGS = -Wf"-dp"

LIBS =  ../libskit.a ../UNSUPP/BLAS1/blas1.o

riters.ex:  riters.o iters.o ilut.o itaux.o $(LIBS)
	$(F77) $(FFLAGS) -o riters.ex riters.o itaux.o $(LIBS)

rilut.ex: rilut.o ilut.o iters.o itaux.o $(LIBS)
	$(F77) $(FFLAGS) -o rilut.ex rilut.o itaux.o $(LIBS)

riter2.ex: riter2.o iters.o ilut.o itaux.o $(LIBS)
	$(F77) $(FFLAGS) -o riter2.ex riter2.o itaux.o $(LIBS)

clean:
	rm -f *.o *.ex core *.trace

../UNSUPP/BLAS1/blas1.o:
	(cd ../UNSUPP/BLAS1; $(F77) $(FFLAGS) -c blas1.f)

../libskit.a:
	(cd ..; $(MAKE) $(MAKEFLAGS) libskit.a)


