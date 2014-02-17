FC = gfortran
FCFLAGS = -g
objects = chemsol.o

main: $(objects)
	$(FC) $(FCFLAGS) -fdump-core -o main main.f90 $(objects)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean

clean: 
	rm -f *.o *.mod *.MOD 
