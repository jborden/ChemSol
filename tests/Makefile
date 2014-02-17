FC = gfortran
FCFLAGS = -g
objects = chemsol.o

main: $(objects)
	$(FC) $(FCFLAGS) -fdump-core -o main main.f90 $(objects)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean test diff

clean: 
	rm -f *.o *.mod *.MOD 

test1: 
	time ./main vdw.par test1.cs >& test1_main.out & 
test2:
	time ./main vdw.par test2.cs >& test2_main.out & 
test3:
	time ./main vdw.par test3.cs >& test3_main.out & 	

diff1: 	
	diff --side-by-side --suppress-common-lines test1.log test1_main.out 
diff2:	
	diff --side-by-side --suppress-common-lines test2.log test2_main.out 
diff3:
	diff --side-by-side --suppress-common-lines test3.log test3_main.out 	
