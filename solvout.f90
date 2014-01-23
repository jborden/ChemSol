subroutine solvout (iterld, do_gas)
  implicit Real*8 (a-h,o-z)
  PARAMETER (MXATM=500)
  parameter (mxcenter=50)
  common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
  common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
  common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
  common /aname/ atom(mxatm),molname,ssname
  character*8 atom
  character*13 molname
  character*4 ssname
  logical do_gas  
  !
  !      elgvn.......noniterative lgvn energy (using distance-dependent dielectric)
  !      elgwa......iterative lgvn energy 

  character*1 dash(72)
  data dash/72*'-'/
  !......................................................................
  write(6,20) dash

  !   -- Total solvation free energy from LD calculation
  erelax = evqdq
  dgsolv=elgwa+etds +evdw+ebw
  if (iterld.eq.0) dgsolv = 0.d0
  dgsolvni=elgvn+ephob+evdw+ebw
  write(6,100)
  if (iterld.eq.0) write(6,85) 
  if (iterld.eq.1) write(6,101) elgwa
  if (iterld.eq.0) write(6,101) elgvn
  write(6,102) etds 
  write(6,103) evdw
  write(6,104) ebw
  if (do_gas) write(6,105) erelax
  !     write(6,107) ecor
  if (iterld.eq.1) write(6,106) dgsolv
  if (iterld.eq.0) write(6,106) dgsolvni
  if (iprint.eq.1) then
     write(6,110)
     write(6,115) elgvn
     write(6,117) dgsolvni
  end if

  !     Short output (do not print for noniterative LD)    

  if (iterld.eq.1) then

     !     Comment the following line for IBM and uncomment the CIBM line
     open (43,file='cs.arc',access='append')
     !IBM  open (43,file='cs.arc',position='append')

     if (do_gas) then
        write(6,199)
        write(6,200) molname,ssname,elgwa, evdw, etds , erelax, ebw, &
             dgsolv-etds,dgsolv
        write(43,201) molname,ssname,elgwa, evdw, etds , erelax, ebw, &
             dgsolv-etds,dgsolv

     else
        write(6,299)
        write(6,300) molname,ssname,elgwa, evdw, etds , ebw, &
             dgsolv-etds,dgsolv
        write(43,301) molname,ssname,elgwa, evdw, etds , ebw, &
             dgsolv-etds,dgsolv
     end if
     close(43)
  end if

199 format(//,'Molecule               lgvn    VdW   -TdS   Relax',&
       '    Born   dHsolv   dGsolv')
200 format (a13,1x,a4,f9.1,3f7.1,3f9.1/)
201 format (a13,1x,a4,f9.1,3f7.1,3f9.1)
299 format(//,'Molecule               lgvn    VdW   -TdS        ',&
       '    Born   dHsolv   dGsolv')
300 format (a13,1x,a4,f9.1,2f7.1,7x,3f9.1/)
301 format (a13,1x,a4,f9.1,2f7.1,7x,3f9.1)

  return
  !......................................................................
20 FORMAT(1X,72A1//' FINAL RESULTS')
85 format(" Noniterative LD results:"/, &
       "Warning! Noniterative LD is not parametrized in cs2.1"/, &
       "         Put iterld=1 in the vdw.par file to run iterative LD"/)
100 format(//,' *****************************'/, &
       ' Langevin Dipoles Calculation '/, &
       ' *****************************'//, &
       ' Transfer of solute at 298K: GAS,1M -> WATER,1M (kcal/mol)',/)
101 format(' Langevin energy         ',10x,f20.2)
102 format(' -TdS                    ',10x,f20.2)
103 format(' VdW energy              ',10x,f20.2)
104 format(' dBorn+dOnsager          ',10x,f20.2)
105 format(' Solute relaxation       ',10x,f20.2)
107 format(' Correlation correction  ',10x,f20.2)
106 format(' Total dG                ',10x,f20.2/)

110 format(//,' ****************************'/,&
       ' Noniterative LD calculation '/,&
       ' ****************************'//,&
       ' Transfer of a solute at 298K: GAS,1M -> WATER,1M',//)
115 format(' Noniterative lgvn energy',10x,f20.2)
117 format(' dG                      ',10x,f20.2)
end subroutine solvout
