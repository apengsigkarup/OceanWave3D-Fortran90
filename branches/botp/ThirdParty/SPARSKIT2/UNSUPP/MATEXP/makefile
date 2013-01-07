FFLAGS = 
F77 = f77

#F77 = cf77
#FFLAGS = -Wf"-dp"

FILES1 = rexp.o exppro.o 
FILES2 = rphi.o phipro.o 

exp.ex:		$(FILES1)
	$(F77) $(FFLAGS) $(FILES1) -o exp.ex  

phi.ex:		$(FILES2)
	$(F77) $(FFLAGS) $(FILES2) -o phi.ex  
