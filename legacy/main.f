! goto 10 50 15 41
C################################################################
C                                                                #
C                        CHEMSOL 2.1                             #
C                                                                #
C               Jan Florian and Arieh Warshel                    #
C                  Department of Chemistry                       #
C             University of Southern California                  #
C                 Los Angeles, CA 90089-1062                     #
C                                                                #
C#################################################################
C                      October 20, 1999     
C#################################################################

C     Program evaluates hydration free energy of neutral and ionic solutes
C     by the Langevin dipole (LD) method.
C     Langevin dipoles are point dipoles positioned at fixed grid points.
C     ChemSol uses MEP point charges to calculate field at grid points.
C     Atomic coordinates are frozen during the calculation.
C     On input: 1/ Cartesian coordinates of the solute /reg1/
C               2/ Electrostatic potential derived atomic charges.
C
C     Changes from the version 1.0:
C  1. To ensure IBM compatibility, initialization of parameters that 
C     appear in common blocks was moved from the SETPAR subroutine to 
C     the BLOCK DATA subprogram.
C  2. Parameter out_cut was decreased from 20 A to 18 A.
C  3. Dimensions of the dipole-dipole pair lists were increased to
C     mxpair=2.5M and mxpair2=5M. About 16MB memory is needed for
C     these dimensions.
C  4. JP3(mxpair) pairlist was eliminated.
C  5. Variable iprint was introduced to reduce printout. By default,
C     iprint=0 (short output).
C  6. A new machine-independent pseudorandom number generator (ran2)
C
C     Changes from the version 1.1
C  1. Bug fix in ran2
C  2. 'implicit double precision' was replaced by 'implicit Real*8'
C     
C     Changes from the version 1.11
C  1. A simple parametrization was done for positive ions from the
C     first and second group of the periodic table. The default 
C     vdW radii have been changed to: Na+(2.55), K+(3.05), Rb+(3.20),
C     Mg++(2.00), Ca++(2.40).
C 
C     Changes from the version 1.12
C  1. New dimensions: mxlgvn=10M, mxpair=5M, mxpair2=10M
C
C     Changes from the version 1.13 
C
C  1. A new 'entropy' function was added (in vlgvn and mu_mu_l) to
C     account for the entropy decrease if solvent dipoles are "frozen"
C     in the regions of the large elstat. field. The sum of the -TdS 
C     and elgvn terms is evaluated to determined the convergence
C     of the iterative process (in sci_lgvn).
C  2. The noniterative LD model is no more supported. (The results
C     are printed only in the main output). This model is used
C     for initiation of the iterative procedure (in lgvnx) as before.
C  3. The Langevin function in vlgvn has been modified to improve
C     the behavior of solvation free energies for the charge separation
C     processes. For details see Table 1S in Suporting Materials
C     for J.Phys.Chem. 1999 paper.
C     The new function increases the contribution of the inner grid and
C     decreases the contribution of the outer grid to the total energy.
C     The factor 3.3 in front of kT represents a screening constant
C     for the field (fj) that is introduced to account for limited
C     dipole-dipole interactions. The new function is the default.
C     The old function can be selected by setting ioldfn=1 in the
C     vdw.par file.
C  4. default vdW radii for some atom types have been changed
C  5. New vdW radii have been added for some transition metals (2+).
C  6. The bulk contribution to the solvation energy is calculated
C     from the molecular dipole and charge (in previous versious it was
C     determined from the total dipole (molecule+LD dipoles))
C     for neutral molecules and from the molecular charge for
C     charged compounds. This modification was introduced to improve
C     treatment of charged compounds with large dipole moments.
C  7. Parameter out_cut was decreased from 18 A to 16 A.
C
C     The version 2.0 of ChemSol program is described in detail
C     in the paper: J. Florian, A. Warshel, "Calculations of Hydration
C     Entropies of Hydrophobic, Polar, and Ionic Solutes in the framework 
C     of the Langevin Dipoles Solvation model",
C     J. Phys. Chem B 1999 (in press). 
C
C     Changes from the version 2.0 

C     The solute relaxation energy is now included implicitly in Elgvn: 
C     dGsolv = dElgvn + dEvdW + dEBorn - TdS 
C     Preffered method for calculation atomic charges is now ESP PCM-B3LYP/6-31G*.


      implicit Real*8 (a-h,o-z)
      PARAMETER (MXATM=500)
      PARAMETER (ONE=1.0d0, ZERO=0.d0)
      parameter (mxcenter=50)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82),
     :             n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)
     :             ,q_gas(mxatm),q_mp2(mxatm)
      common /aname/ atom(mxatm),molname,ssname
      common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
      common /born/ rg_reg1, rgim
      common /pcgrid/ drg,rg,dxp0(3),
     :              rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /pctimes/ ndxp,itl,itp
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /lra/ clgvn, slgvn
      character*8 atom,dumm1 
      character*8 molname
      character*4 ssname
      character*256 fname
      character*256 input_file_name
      character*256 vdw_name
      character*3 relax 
      logical do_gas
      integer iac_conv (89), iacp(mxatm), mass(88)
      data iac_conv/1,2,
     1              3,4,5,6,9,13,16,17,
     2              18,19,20,21,22,24,25,26,
     3     27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,
     4     45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,
     5     63,64,65,
     6     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     7              66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,
     8     81,82/
c      ChemSol atom types: 
c       1-H0, 2-He, 3-Li, 4-Be, 5-B, 6-C0, 7-C1, 8-C2, 9-N0, 10-N1,
c       11-N2, 12-X, 13-O0, 14-O1, 15-O2, 16-F, 17-Ne, 18-Na, 19-Mg, 20-Al,
c       21-Si, 22-P, 23-X, 24-S, 25-Cl, 26-Ar, 27-K, 28-Ca, 29-Sc, 30-Ti, 
c       31-V, 32-Cr, 33-Mn, 34-Fe, 35-Co, 36-Ni, 37-Cu, 38-Zn, 39-Ga,40-Ge,
c       41-As, 42-Se,43-Br, 44-Kr, 45-Rb, 46-Sr, 47-Y, 48-Zr, 49-Nb, 50-Mo, 
c       51-Tc, 52-Ru, 53-Rh, 54-Pd, 55-Ag, 56-Cd, 57-In, 58-Sn, 59-Sb, 60-Te, 
c       61-I, 62-Xe, 63-Cs, 64-Ba, 65-La, 66-Hf, 67-Ta, 68-W, 69-Re, 70-Os, 
c       71-Ir, 72-Pt, 73-Au, 74-Hg, 75-Tl, 76-Pb, 77-Bi, 78-Po, 79-At, 80-Rn, 
c       81-Fr, 82-Ra.

      data mass /1, 4,
     1           7,9,11,12,14,16,19,20,
     2           23,24,27,28,31,32,35,40,
     3           39,40,45,48,51,52,55,56,59,59,64,65,70,73,75,79,80,84,
     4 85,88,89,91,93,96,97,101,103,106,108,112,115,119,122,128,127,131,
     5 132,137,139,
     6 140,141,144,145,150,152,157,159,163,165,167,169,173,175,
     7 178,181,184,186,190,192,195,197,201,204,207,208,209,210,222,
     8 223,226/ 

c     open archive (this is actually done in subroutine solvout
c     open (43,file='cs.arc',access ='append')
c     open input file with vdw and grid options
c     open (44, file='vdw.par')
      call getarg(1,vdw_name)
      open (44,file=vdw_name)
c     open input file for atom input
c     call getenv('SOLVINP',fname)
      call getarg(2,input_file_name)
      open (45,file=input_file_name)
      read (45,'(a13)') molname
      read (45,*) n_reg1, ngeom
      imp2 = 0 
      read(45,'(a4)') ssname
      do i=1,n_reg1
      read (45,1000) atom(i),zan(i),q(i),xw(1,i), xw(2,i), xw(3,i)
c     Gaussian atom types (nuclear charge) are mapped onto ChemSol ones (iacw). 
      iacp(i)=int(zan(i))
      iacw(i)=iac_conv(iacp(i))
      latom(i)=i
      end do

c     Calculate the molecular mass (a.u.)
      amas = 0.0d0
      do i=1,n_reg1
      amas =amas + dble(mass(int(zan(i))))
      end do

c     If charges from a PCM calculation (q_pcm) are available, they will
c     be used in the evaluation of dGlgvn. Charges can be inputed in two ways:
c     1/ q_pcm are given as the first set of charges in the
c     input file and no second set of charges is given. This scheme
c     is useful when calculating charges using the Gaussian98 code
c     to avoid extra calculation of the gas-phase charges.
c     2/ q_gas charges are given as the first set of charges and q_pcm
c     are given as the second set of charges. This is the input format 
c     used in Chemsol 1.0 - 2.0. The advantage of this scheme is that
c     the additional information about the gas-phase charge distribution 
c     can be used to evaluate EXPLICITLY the contribution of the solute polarization
c     by the solvent to the total solvation free energy (dGsolv). This 
c     contribution is called here dGrelax. Because dGrelax is implicitly
c     included in dGsolv if pcm charges are available, dGsolv results
c     are not affected by the presence/absence of the gas_phase charges 
C     in the input file.


      do_gas = .false.
      read (45,'(a3)') relax
      if (relax.eq.'pcm') then
        do_gas = .true.
        do i=1,n_reg1
          q_gas(i) = q(i)
        end do

        do i=1,n_reg1
          read (45,'(17x,f10.4)') q(i)
        end do
      else
        backspace(45)
        do i=1,n_reg1
          q_gas(i) = q(i)
        end do
      end if

c     If charges from the correlated calculation are available, they will
c     be used in the evaluation of dGmp2.
c     Note that in cs2.1 imp2 is always zero.
c     The unused code is still kept in case a construct like this is needed in
c     future.
      if(imp2.eq.1) then
      read (45,'(a3)') corr
      do i=1,n_reg1
      read (45,'(17x,f10.4)') q_mp2(i)
      end do
      else
      do i=1,n_reg1
      q_mp2(i) = q(i)
      end do
      end if
      
c     Check consistency of charges
      qsum1 = 0.
      qsum2 = 0.
      qsum3 = 0.
      do i=1,n_reg1
      qsum1 = qsum1 + q(i)
      qsum2 = qsum2 + q_gas(i)
      qsum3 = qsum3 + q_mp2(i)
      end do

c     Calculate translational entropy of the gas-phase solute (1M ideal
c     gas, T=298.5)
      weight = amas
      TdS_gas = entropy (weight)   
      
 55   format(
     $6x,"***********************************************************"/
     $6x,"                                                           "/
     $6x,"                       CHEMSOL 2.1                         "/
     $6x,"                                                           "/
     $6x,"              Jan Florian and Arieh Warshel                "/
     $6x,"            University of Southern California              "/
     $6x,"                    Los Angeles, 1999                      "/
     $6x,"                                                           "/
     $6x,"***********************************************************"//
     $      )

      write (6,55)
      write (6,'(/" Solute                           : ",
     $           5x,a13)')molname
      write (6,'(" Total charge for gas-phase solute: ",f9.2)') qsum1
      write (6,'(" Total charge for PCM solute      : ",f9.2)') qsum2
c     write (6,'(" Total MP2 charge                 : ",f9.2)') qsum3
      write (6,'(" Molecular weight (a.u.)          : ",f9.2)') amas
      write (6,'(" Gas-phase transl. entropy (298dS): ",f9.2,      
     $            " kcal/mol")') TdS_gas
      k1 = 0
      k2 = 0
      k3 = 0
      do i=1,10
      if (abs( abs (qsum1)-float(i)+1.0).lt. 0.012) k1=1
      if (abs (abs (qsum2)-float(i)+1.0).lt. 0.012) k2=1
      if (abs (abs (qsum3)-float(i)+1.0).lt. 0.012) k3=1
      end do
      if (k1.eq.0) stop "Inconsistent gas-phase charges"
      if (k2.eq.0) stop "Inconsistent pcm charges"
      if (k3.eq.0) stop "Inconsistent correlated charges"
      
c      The second character of the atom name is used to make wider selection 
c      of ChemSol atom types (that differ in VdW radii).
      do i=1,n_reg1 
      if (iacw(i).ne.1) then
        if (atom(i)(2:2) .eq. '1') iacw(i) = iacw(i)+1
        if (atom(i)(2:2) .eq. '2') iacw(i) = iacw(i)+2
        if (atom(i)(2:2) .eq. '3') iacw(i) = iacw(i)+3
      end if
      if (iacw(i) .eq.0) then
      write(6,'("Unknown atom type on input:",i3,a8,f4.1)') i,atom(i),zan(i)
      stop 
      end if
      end do

      do 30 j=1,ngeom
      if (j.ge.2) then
         write(6,2000) j
         read(45,'(a4)') ssname
         do i=1,n_reg1
         read (45,1000) dumm1,dumm2,q(i),xw(1,i), xw(2,i), xw(3,i)
         end do
         read (45,'(a3)') relax
         if (relax.eq.'pcm') then
            do i=1,n_reg1
            read (45,'(17x,f10.4)') q_gas(i)
            end do
            else
            do i=1,n_reg1
            q_gas(i) = q(i)
            end do
         end if
         if (imp2.eq.1) then
            read (45,'(a3)') corr
            do i=1,n_reg1
            read (45,'(17x,f10.4)') q_mp2(i)
            end do
            else
            do i=1,n_reg1
            q_mp2(i) = q(i)
            end do
         end if
         open (44, file='vdw.par')
      end if

c     Initialize parameters and read in the option file (vdw.par).
      call readopt (iterld)

c     Calculate solvation (LD method).
c     Set parameter iprint to 1, if more output is needed (for debug)
      iprint = 0
      call dg_ld (iterld,iprint)

c     write out solvation energy
      call solvout (iterld,do_gas)

      close (44)
30    continue

      close (43)
      close (45)
1000  format(1x,A8,F8.1,2F10.4,3F9.4)
2000  format(i3)
      end
