FC = gfortran
FCFLAGS = -g 
objects = chemsol.o push_array.o


cs: $(objects)
	$(FC) $(FCFLAGS) -fdump-core -o chemsol main.f90 $(objects) 	

cs21: 
	$(FC) $(FCFLAGS) -o cs21 legacy/cs21.f

chemsol.o: push_array.o

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean

clean: 
	rm -f *.o *.mod *.MOD 
