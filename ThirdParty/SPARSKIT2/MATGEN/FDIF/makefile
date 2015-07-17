FFLAGS = 
F77 = f77

#F77 = cf77
#FFLAGS = -Wf"-dp"

FILES1 = rgen5pt.o functns.o 

gen5.ex:	$(FILES1) ../../libskit.a 
	$(F77) $(FFLAGS) -o gen5.ex $(FILES1) ../../libskit.a 

FILES2 = rgenblk.o functns.o 
genbl.ex:	$(FILES2)  ../../libskit.a 
	$(F77) $(FFLAGS) -o genbl.ex $(FILES2) ../../libskit.a 

clean:
	rm -f *.o *.ex core *.trace

../../libskit.a:
	(cd ../..; $(MAKE) $(MAKEFLAGS) libskit.a)

