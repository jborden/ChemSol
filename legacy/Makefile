FC = gfortran
FCFLAGS = -g

cs21:	
	$(FC) $(FCFLAGS) -fdump-core -o cs21 cs21.f

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<
%.o: %.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean test diff

clean: 
	rm -f *.o *.mod *.MOD 
