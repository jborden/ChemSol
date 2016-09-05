# ChemSol 
## Arieh Warshel, Jan Florian,  James Borden

ChemSol is a program for the calculation of solvation energies by using the Langevin Dipole model of the solvent and ab-initio calculations. The advantages of it over other solvation models is that it can more accurately model solvation effects for molecules with a high net charge, i.e. over +1. It has been particuarly parametrized for reactions involving organic reactions with phosphorous containing compounds such as organophosphates by water.

However, it utilizes the Langevin Dipole model that is a generic solvation model that is flexible enough to model non-aqueous solutions and reactions involving inorganic compounds.

# Modification of ChemSol 2.1

The original code for Chemsol, hosted at http://laetro.usc.edu/software.html, is limited in the size of the system it can model. The size of the system is defined as the total amount of atoms. It has been modified by James Borden to model the solvation free energy of much larger systems such as proteins and enzymes. It has also been updated to use F9X/200X instead of F77 which it was original coded in. The code has been re-factored to use Fortran modules and programmed in a functional style. 

# Compilation

A makefile has been created that utilizes the GNU gfortran compiler.  If you would like to compile it using ifort, simply modify the FC variable from 'gfortran' to the compiler you wish to use. e.g. 'ifort'

To compile the code, execute

```bash
$ make chemsol
```

The original code in the legacy/ dir can be compiled by typing

```bash
$ make legacy
```
If you wish to use compile a version that is debuggable with gdb, type

```bash
$ make debug
```

# Running ChemSol

Documentation on for setting up calculations and parametrization can be found in documentation/cs21_manual.pdf. The new version of chemsol can be run with the folliwing command:

```bash
$ chemsol vdw.par model.cs >& model.out &
```

# Tests

The test dir contains output from legacy chemsol 2.1 to compare with the new version. A test can be run with the command 

```bash
$ make test1
```

and compared to the old output with

```bash
$ make diff1
```

Please see tests/Makefile for more examples

# Large Systems

Using ChemSol for larger systems has not been fully tested. If you would like to use it to calculate the solvation free energies of larger systems, please contact me (James Borden) via GitHub. I will gladly work with you and Jan Florian to make sure that your system has correct solvation free energies calculated for it. 
# Questions and Concerns

Please do not hesitate to contact the authors if you need help or assitance. Feel free to open an issue or to make pull requests!



