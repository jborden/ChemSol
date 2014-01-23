subroutine ran_shift(i,center1,center2)
  implicit Real*8 (a-h,o-z)

  parameter (mxcenter=50)
  PARAMETER (MXATM=500)
  common /pctimes/ ndxp,itl,itp
  common /pcgrid/ drg,rg,dxp0(3),&
       rg_inner,drg_inner
  common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
  dimension oshift(3*mxcenter)
  dimension dxp(3), center2(3), center1(3)
  save oshift

  ! --   Initialize the random number generator and
  !      generate random origin shifts for ndxp grids.
  if (i.eq.1) then
     iseed = -931
     idum = 1
     dumm = ran2(iseed)
     do kk = 1, 3*ndxp
        oshift(kk) = ran2 (idum)
     end do
  end if

  fact=drg_inner
  if(rg_inner.eq.0.d0) fact=drg
  if(i.eq.1)fact=0.0d0
  dxp(1)=fact*(1.d0-2.d0*oshift(3*i-2))
  dxp(2)=fact*(1.d0-2.d0*oshift(3*i-1))
  dxp(3)=fact*(1.d0-2.d0*oshift(3*i))

  if(i.ne.1)then
     center2(1)=center1(1)+dxp(1)+dxp0(1)
     center2(2)=center1(2)+dxp(2)+dxp0(2)
     center2(3)=center1(3)+dxp(3)+dxp0(3)
  else
     center2(1)=center1(1)+dxp0(1)
     center2(2)=center1(2)+dxp0(2)
     center2(3)=center1(3)+dxp0(3)
  endif

  temp_center(i,1)=center2(1)
  temp_center(i,2)=center2(2)
  temp_center(i,3)=center2(3)

  !      write(6,100) center1, dxp0, dxp, center2
  write(6,100) center2

  return
  !......................................................................
100 format(/ &
         !    s ' original grid origin  ',3f9.3/
         !    s ' original origin shift ',3f9.3/
         !    s ' random origin shift   ',3f9.3/
       ' Grid origin            ',3f9.3)
end subroutine ran_shift
