program main
  use chemsol, only: entropy,readopt,dg_ld,solvout 
  integer,parameter  :: mxlgvn = 10000 ! maximum allowed langevin dipoles,  should use dynamic arrays in gen_gridx
  integer,parameter :: mxcenter = 50 ! needed by ran_shift and elgvn_ave
  !      ChemSol atom types: 
  !       1-H0, 2-He, 3-Li, 4-Be, 5-B, 6-C0, 7-C1, 8-C2, 9-N0, 10-N1,
  !       11-N2, 12-X, 13-O0, 14-O1, 15-O2, 16-F, 17-Ne, 18-Na, 19-Mg, 20-Al,
  !       21-Si, 22-P, 23-X, 24-S, 25-Cl, 26-Ar, 27-K, 28-Ca, 29-Sc, 30-Ti, 
  !       31-V, 32-Cr, 33-Mn, 34-Fe, 35-Co, 36-Ni, 37-Cu, 38-Zn, 39-Ga,40-Ge,
  !       41-As, 42-Se,43-Br, 44-Kr, 45-Rb, 46-Sr, 47-Y, 48-Zr, 49-Nb, 50-Mo, 
  !       51-Tc, 52-Ru, 53-Rh, 54-Pd, 55-Ag, 56-Cd, 57-In, 58-Sn, 59-Sb, 60-Te, 
  !       61-I, 62-Xe, 63-Cs, 64-Ba, 65-La, 66-Hf, 67-Ta, 68-W, 69-Re, 70-Os, 
  !       71-Ir, 72-Pt, 73-Au, 74-Hg, 75-Tl, 76-Pb, 77-Bi, 78-Po, 79-At, 80-Rn, 
  !       81-Fr, 82-Ra.
  integer,dimension(89),parameter :: iac_conv = [1,2, &
       3,4,5,6,9,13,16,17, &
       18,19,20,21,22,24,25,26, &
       27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44, &
       45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62, &
       63,64,65, &
       0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
       66,67,68,69,70,71,72,73,74,75,76,77,78,79,80, &
       81,82]

  integer,dimension(88),parameter :: mass = [1,4, &
       7,9,11,12,14,16,19,20, &
       23,24,27,28,31,32,35,40, &
       39,40,45,48,51,52,55,56,59,59,64,65,70,73,75,79,80,84, &
       85,88,89,91,93,96,97,101,103,106,108,112,115,119,122,128,127,131, &
       132,137,139, &
       140,141,144,145,150,152,157,159,163,165,167,169,173,175, &
       178,181,184,186,190,192,195,197,201,204,207,208,209,210,222, &
       223,226]

  ! core topology
  character(3),allocatable :: atom(:)  ! Chemical Symbol
  real(8),allocatable :: zan(:),q(:),xw(:,:) ! atomic number, atomic charge, cartesian coordinates of nucleus
  ! Van der Waal radii topology
  real(8),allocatable :: rpi(:)
  ! optional charge topology
  real(8),allocatable :: q_gas(:),q_mp2(:)
  ! chemsol atom type topology
  integer,allocatable :: iacw(:)

  real(8),dimension(82) :: &
       !      H   He  Li   Be   B    C    C1  C2   N    N1
       rp = [2.3,2.3,2.15,2.00,2.40,2.65,3.0,3.25,2.65,2.85, &
       !N2   X   O   O1  O2    F    Ne  Na   Mg  Al
       3.2,3.0,2.32,2.65,2.8,2.46,2.5,2.58,1.82,1.70, &
       !Si  P    X  S    Cl  Ar   K    Ca   Sc  Ti
       3.1,3.2,3.0,3.2,3.16,2.8,3.06,2.38,1.5,2.00, &
       !V    Cr   Mn   Fe   Co   Ni  Cu1+  Zn  Ga   Ge
       1.91,1.89,1.93,1.84,1.57,1.50,1.88,1.55,2.00,2.50, &
       !As   Se    Br   Kr   Rb  Sr    Y    Zr  Nb    Mo
       3.00,3.00,3.44,3.00,3.25,2.70,1.75,2.00,2.00,2.00, &
       !Tc    Ru   Rh  Pd   Ag    Cd   In   Sn   Sb  Te
       2.00,2.00,2.00,2.00,2.25,1.98,2.45,2.44,3.00,3.70, &
       !I    Xe   Cs   Ba   La   Hf   Ta   W    Re   Ir 
       3.80,3.55,3.58,2.92,2.30,3.00,3.00,3.00,3.00,3.00, &
       !Pt   Au   Hg   Tl   Pb   Bi   Po   At   Rn   Fr
       3.00,1.77,1.94,2.00,2.55,3.00,3.00,3.00,3.00,3.00, &
       !Ra   Ac
       3.50,3.00]

  real(8) :: vdwc6(82)
  ! ecor was removed as it seems to be local to dg_ld
  real(8) :: ebw,elgvn,elgwa,erelax,evdw,evqdq
  real(8) :: rg_reg1, rgim
  real(8) :: drg,rg,dxp0(3),rg_inner,drg_inner
  real(8) :: rdcutl,out_cut
  integer :: ndxp,itl,itp
  real(8) :: ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
  real(8) :: vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
  real(8) :: clgvn, slgvn
  real(8) :: pcenter(3)


  character(8):: dumm1
  character(13) :: molname
  character(4) :: ssname
  character(256) :: fname
  character(256) :: input_file_name
  character(256) :: vdw_name
  character(3) :: relax 
  logical :: do_gas
  !integer ::iac_conv (89)
  integer :: n_reg1,n_inner
  ! from blockdata
  rg = 26.0
  drg = 3.0
  drg_inner = 1.0
  rdcutl = 6.1
  out_cut = 16.0
  dxp0 = [0.0,0.0,0.0]
  rgim = 2.0
  !ndxp = 10
  itl = 399
  itp = 5
  rzcut = 1.1
  !     open archive (this is actually done in subroutine solvout
  !     open (43,file='cs.arc',access ='append')
  !     open input file with vdw and grid options
  !     open (44, file='vdw.par')
  call getarg(1,vdw_name)
  !write (*,*) trim(vdw_name)
  open (44,file=vdw_name)
  !     open input file for atom input
  !     call getenv('SOLVINP',fname)
  call getarg(2,input_file_name)
  !      write (*,*) input_file_name
  open (45,file=input_file_name)
  read (45,'(a13)') molname
  read (45,*) n_reg1, ngeom
  imp2 = 0 
  read(45,'(a4)') ssname
  
  allocate(atom(n_reg1),zan(n_reg1),q(n_reg1),xw(3,n_reg1))
  allocate(iacw(n_reg1))
  do i=1,n_reg1
     ! read (45,1000) atom(i),zan(i),q(i),xw(1,i), xw(2,i), xw(3,i)
     ! just read in the line, no need for special formatting
     read (45,*) atom(i),zan(i),q(i),xw(1,i), xw(2,i), xw(3,i)
     !     Gaussian atom types (nuclear charge) are mapped onto ChemSol ones (iacw). 
     ! write (*,*) atom(i),zan(i),q(i),xw(1,i),xw(2,i),xw(3,i)
     iacw(i)=iac_conv(int(zan(i)))
  end do
  !     Calculate the molecular mass (a.u.)
  amas = 0.0d0
  do i=1,n_reg1
     amas = amas + dble(mass(int(zan(i))))
  end do

  !     If charges from a PCM calculation (q_pcm) are available, they will
  !     be used in the evaluation of dGlgvn. Charges can be inputed in two ways:
  !     1/ q_pcm are given as the first set of charges in the
  !     input file and no second set of charges is given. This scheme
  !     is useful when calculating charges using the Gaussian98 code
  !     to avoid extra calculation of the gas-phase charges.
  !     2/ q_gas charges are given as the first set of charges and q_pcm
  !     are given as the second set of charges. This is the input format 
  !     used in Chemsol 1.0 - 2.0. The advantage of this scheme is that
  !     the additional information about the gas-phase charge distribution 
  !     can be used to evaluate EXPLICITLY the contribution of the solute polarization
  !     by the solvent to the total solvation free energy (dGsolv). This 
  !     contribution is called here dGrelax. Because dGrelax is implicitly
  !     included in dGsolv if pcm charges are available, dGsolv results
  !     are not affected by the presence/absence of the gas_phase charges 
  !     in the input file.

  allocate(q_gas(n_reg1))

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

  !     If charges from the correlated calculation are available, they will
  !     be used in the evaluation of dGmp2.
  !     Note that in cs2.1 imp2 is always zero.
  !     The unused code is still kept in case a construct like this is needed in
  !     future.
  allocate(q_mp2(n_reg1))

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

  !     Check consistency of charges
  qsum1 = 0.
  qsum2 = 0.
  qsum3 = 0.
  do i=1,n_reg1
     qsum1 = qsum1 + q(i)
     qsum2 = qsum2 + q_gas(i)
     qsum3 = qsum3 + q_mp2(i)
  end do

  !     Calculate translational entropy of the gas-phase solute (1M ideal
  !     gas, T=298.5)
  !      weight = amas
  !      TdS_gas = entropy(amas)   

55 format(6x,"***********************************************************"/ &
       6x,"                                                           "/ &
       6x,"                       CHEMSOL 2.1                         "/ &
       6x,"                                                           "/ &
       6x,"              Jan Florian and Arieh Warshel                "/ &
       6x,"            University of Southern California              "/ &
       6x,"                    Los Angeles, 1999                      "/ &
       6x,"                                                           "/ &
       6x "***********************************************************"/)

  write (6,55)
  write (*,*)
  write (6,'(/" Solute                           : ",5x,a13)') molname
  write (6,'(" Total charge for gas-phase solute: ",f9.2)') qsum1
  write (6,'(" Total charge for PCM solute      : ",f9.2)') qsum2
  !     write (6,'(" Total MP2 charge                 : ",f9.2)') qsum3
  write (6,'(" Molecular weight (a.u.)          : ",f9.2)') amas
  write (6,'(" Gas-phase transl. entropy (298dS): ",f9.2," kcal/mol")') entropy(amas)
  k1 = 0
  k2 = 0
  k3 = 0
  ! essentially checks to make sure that the charges are whole numbers within the range
  ! of 9 to -9
  do i=1,10
     if (abs (abs (qsum1)-float(i)+1.0).lt. 0.012) k1=1
     if (abs (abs (qsum2)-float(i)+1.0).lt. 0.012) k2=1
     if (abs (abs (qsum3)-float(i)+1.0).lt. 0.012) k3=1
  end do
  if (k1.eq.0) stop "Inconsistent gas-phase charges"
  if (k2.eq.0) stop "Inconsistent pcm charges"
  if (k3.eq.0) stop "Inconsistent correlated charges"

  !      The second character of the atom name is used to make wider selection 
  !      of ChemSol atom types (that differ in VdW radii).
  do i=1,n_reg1 
     if (iacw(i).ne.1) then
        ! write (*,*) atom(i)(2:2) debug
        if (atom(i)(2:2) .eq. '1') iacw(i) = iacw(i)+1
        if (atom(i)(2:2) .eq. '2') iacw(i) = iacw(i)+2
        if (atom(i)(2:2) .eq. '3') iacw(i) = iacw(i)+3
     end if
     if (iacw(i) .eq.0) then
        write(6,'("Unknown atom type on input:",i3,a8,f4.1)') i,atom(i),zan(i)
        stop 
     end if
  end do

  !      do 30 j=1,ngeom
  do j=1,ngeom
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

     !     Initialize parameters and read in the option file (vdw.par).
     !call readopt (iterld)
     iprint = 0
     allocate(rpi(n_reg1))
     call readopt(iterld,vdwc6,dxp0,clgvn,slgvn,tds0,rp,vdwsl,phobsl,ephil1,ephil2,rzcut, &
          rpi,pcenter,rg_reg1,rg,rg_inner,rgim,ndxp,iacw,xw,q,q_gas,n_reg1,drg,iprint)

     !     Calculate solvation (LD method).
     !     Set parameter iprint to 1, if more output is needed (for debug)
     !     call dg_ld (iterld,iprint)
     call dg_ld (iterld,iprint,evqdq,elgwa,evdw,elgvn,ephob,etds,ebw,clgvn,dxp0,ephil1,ephil2,iacw,ndxp, & 
          pcenter,phobsl,q,q_gas,q_mp2,rg,rg_inner,rg_reg1,rgim,rpi, &
          rzcut,slgvn,tds0,vdwc6,vdwsl,xw,atom,n_reg1, &
          drg,drg_inner,rdcutl,out_cut,itl)

     !     write out solvation energy
     !call solvout (iterld,do_gas)
     call solvout (iterld,iprint,do_gas,evqdq,elgwa,etds,evdw,ebw,elgvn,ephob,molname,ssname)
     close (44)
     ! deallocate molecular topology
     deallocate(atom,zan,q,xw,rpi,q_gas,q_mp2,iacw)
     !30    continue
  end do
  close (43)
  close (45)
1000 format(1x,A3,F8.1,2F10.4,3F9.4)
2000 format(i3)
end program main
