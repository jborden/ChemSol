FC = gfortran
FCFLAGS = -g -O2
objects = BLOCKDATA1.o dg_ld.o ef_ld.o elgvn_ave.o entropy.o gen_gridx.o lgvnx.o mu_mu_l.o newf_lcut.o pairlistw.o ran2.o ran_shift.o readopt.o sci_lgvn.o solvout.o updatelong.o vatom.o vbornx.o vlgvn.o

main: $(objects)
	$(FC) $(FCFLAGS) -o main main.f $(objects)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean

clean: 
	rm -f *.o *.mod *.MOD
