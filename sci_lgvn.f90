subroutine sci_lgvn(energy,elgvni,tds,nd,icent)
! --  Iterative calculation of langevin dipoles.
      implicit Real*8 (a-h,o-z)
      parameter (mxatm=500)
      parameter (mxlgvn=10000)
!     parameter (mxpair3=1800000)
      parameter (mxcenter=50)
      parameter (mxpair2=10000000)
      parameter (mxpair=5000000)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
           n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm), &
           q_gas(mxatm),q_mp2(mxatm)
      common /pcefaefb/efa(3,mxlgvn)
      common efal (3,mxlgvn)
      common /pcdipcut/ rdcutl,out_cut
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn), &
           da(3,mxlgvn)
      common /pctimes/ ndxp,itl,itp
      integer*2 jp, jp2
      common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn), &
           jp2(mxpair2),ip3(0:mxlgvn)
      common /lra/ clgvn, slgvn

!:::  input vars
      real*8 energy
      integer nd
 
!:::  local vars
      integer i,l,icent,iopen
      real*8 conv2
      real*8 epom(44)
!     integer*2 jp3(mxpair3)
      character*1 dash(72)
      data dash/72*'-'/
      data iopen/1/
!......................................................................
      tds = 0.0d0
      write(6,1006) 
!     write(6,1002) nd
!     call pairlistw(nd,xd,jp3)
      call pairlistw(nd,xd)
      write(6,1008)
      conv2=-332.d0*slgvn
      do i=1,10
      epom(i) = 0.d0
      end do
! --  Relaxation of the langevin dipoles
      do 40 l=1,itl

!     Dump dipoles into xyz file for viewing by the program xmol
      if (.false.) then
      open (50, file='cs_dipoles.xyz')
      nd2plot=0
      do i=1,nd
         if (xl(1,i).eq.0.0) nd2plot=nd2plot+1
      enddo
      write(50,2006) nd2plot
      write(50,*)
      do i=1,nd
         scale = 25.
         if (i.le.n_inner) scale = 50.
         if (xl(1,i).eq.0.0) then
         write (50,2007) xl(1,i),xl(2,i), xl(3,i), &
              scale*xmua(1,i),scale*xmua(2,i),scale*xmua(3,i)
         endif
      enddo
      end if
 2006 format(i5)
 2007 format('X',6f10.3)


      eita=0.
      do i=1,n_inner
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eiti=conv2*eita
      do i=n_inner+1,nd
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eita=conv2*eita
      write(6,1013) l-1,eita,eiti,eita-tds
!     Remember 0th iteration energy (noniterative lgvn energy).
      if (l.eq.1) elgvni=eita

! --  Calculate deviation of the last computed free energy from the average
!     of previous ni iterations. If this is less than 0.1% langevin
!     energy is converged. 
      ni=10
      if(l.gt.5) then
        eave=0.d0
        epom(ni+1)= eita - tds
        do i=1,ni
        epom(i) = epom(i+1)
        eave=eave+epom(i)
        end do

        eave=eave/float(ni)
        iconvp=1
        iconv=1
        do i=1,ni
        if((100.d0*abs((eave-epom(i))/eave)).gt.0.1d0) iconv=0
        iconvp = iconvp * iconv
        end do

! --  Alternatively, stop iterations if the minimum of the dG is reached.
        if(l.gt.16.and.epom(ni).gt.epom(ni-1).and. &
             epom(ni-1).gt.epom(ni-2))  iconvp = 1
        if(iconvp.eq.1) then 
!       write(6,'("Elgvn converged. Average = ",f8.3)') eave
! #GOTO 50
        goto 50
        end if
      end if

      call newf_lcut(nd,l)
! -- In the local reaction field approximation, long-range dipole-dipole
!    interactions are updated only in the first n steps; the field from
!    distant dipoles is fixed at its last value afterwards.
!     if(l.lt.10.or.l.eq.45) call updatelong(nd,l,jp3)
      if(l.lt.10.or.l.eq.45) call updatelong(nd,l)
      stepaw=0.2
      if(l.lt.4) stepaw=0.10
      if(l.lt.11) stepaw=0.20
      if(l.gt.28) stepaw=0.30
      if(l.gt.80) stepaw=0.50
      call mu_mu_l (nd,stepaw,tds)
 40   continue

      eita=0.d0
      do i=1,n_inner
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eiti=conv2*eita
      do i=n_inner+1,nd
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eita=conv2*eita
      write(6,1014) l-1,eita,eiti,eita-tds

 50   continue
      energy=eita

      return
!....................................................................
 1002 format(1x,'Total number of dipoles = ',i8)
 1006 format(/,1x,'Iterative LD calculation for the chosen grid')
 1008 format(//1x,'Iteration',9x,'Elgvn',5x,'Elgvn_inner',5x, &
           'Elgvn-TdS_immob'/)
 1013 format(5x,i3,4x,f12.3,5x,f12.3,5x,f12.3)
 1014 format(5x,i3,4x,f12.3,5x,f12.3,5x,f12.3//)
    end subroutine sci_lgvn
