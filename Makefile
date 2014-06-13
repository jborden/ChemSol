FC = gfortran
FCFLAGS = -g 
objects = chemsol.o push_array.o

main: $(objects)
	$(FC) $(FCFLAGS) -fdump-core -o main main.f90 $(objects) 

cs21: 
	$(FC) $(FCFLAGS) -o cs21 cs21.f

chemsol.o: push_array.o

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean

clean: 
	rm -f *.o *.mod *.MOD 
