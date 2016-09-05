FC = gfortran
FCFLAGS =
SRCDIR = src
objects = $(SRCDIR)/chemsol.o $(SRCDIR)/push_array.o

chemsol: $(objects)
	$(FC) $(FCFLAGS) -o chemsol $(SRCDIR)/main.f90 $(objects) 	

debug: $(objects)
	$(FC) -g -o chemsol $(SRCDIR)/main.f90 $(objects)

legacy: 
	$(FC) $(FCFLAGS) -o cs21 legacy/cs21.f


$(SRCDIR)/chemsol.o: $(SRCDIR)/push_array.o

$(SRCDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(SRCDIR)/%.o: %.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean legacy

clean: 
	rm -rf *.o *.mod *.MOD chemsol cs21 chemsol.dSYM $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.MOD
