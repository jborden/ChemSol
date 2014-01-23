subroutine elgvn_ave(iterld,ncenter)
  implicit Real*8 (a-h,o-z)
  parameter (mxcenter=50)
  common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
  common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
  common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
  common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
  common /volume/ nvol(mxcenter)

  character*1 dash(72)
  data dash/72*'-'/
  !...................................................................
  eav=0.d0
  eboltz=0.d0
  enorm=0.d0
  write(6,1000) 
  temper=300.d0

  wlowest=10000.d0
  do  jj=1,ncenter
     if(temp_elgvn(jj).lt.wlowest) then
        wlowest=temp_elgvn(jj)
        nlowest=jj
     endif
  end do

  elowest=temp_elgvn(nlowest)
  do  i=1,ncenter
     write(6,1001) i,temp_elgvn(i)
     eav=eav+temp_elgvn(i)/dble(ncenter)
     expe=dexp(-(temp_elgvn(i)-elowest)/0.002d0/temper)
     eboltz=eboltz+temp_elgvn(i)*expe
     enorm=enorm+expe
  end do

  eboltz=eboltz/enorm
  write(6,1002)elowest,eav,eboltz,temper

  sum_vdwl =0.0
  sum_phob =0.0
  sum_tds  =0.0
  sum_surf =0.0
  sum_vol = 0.0
  do j=1,ncenter
     sum_vdwl = sum_vdwl + evdwl(j)
     sum_surf = sum_surf + fsurfa(j)
     sum_phob = sum_phob+phobsl*fsurfa(j)
     sum_tds = sum_tds + tdsl(j)
     sum_vol = sum_vol + dble(nvol(j))
  enddo
  ave_surf = sum_surf/dble(ncenter)
  evdw = sum_vdwl/dble(ncenter)
  ephob = sum_phob/dble(ncenter)
  etds = sum_tds/dble(ncenter)
  vol = sum_vol/dble(ncenter)


  ! --  Add solute free-volume entropy and hydrophobic entropy:
  write(6,'(/," Contributions to -TdS (kcal/mol)")' )
  write(6,'('' Change in free-volume:  '',5x,f10.3)')  tds0
  write(6,'('' Hydrophobic          :  '',5x,f10.3)')  ephob
  write(6,'('' Dipolar saturation   :  '',5x,f10.3)')  -etds

  !     write(6,'(/,''Solute volume (A**3):  '',f10.3)') vol 

  etds = tds0 + ephob - etds

  if (iprint.eq.1) then
     write(6,'(/,''Average VDW energy:            '',f10.3)') evdw
     write(6,'(''Average -TdS (kcal/mol):       '',f10.3)') etds
     write(6,'(''Average solvation entropy (eu):'',f10.3)') &
          -1000.d0*etds/temper
     write(6,'(''Average relaxation energy:     '',f10.3)') evqdq
  end if

  if (iterld.eq.1) elgwa = eav 
  !....................................................................
1000 format(1x,' Final LD energies for different grids'// &
       '                    Elgvn     '/)
1001 format(1x,i10,2x,f10.3)
1002 format(//' lowest energy        --->  ',f10.3/ &
       ' mean energy          --->  ',f10.3/ &
       ' boltzmann <energies> --->  ',f10.3/ &
       ' for temperature      --->  ',f10.1,' K'/)
  return
end subroutine elgvn_ave
