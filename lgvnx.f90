subroutine lgvnx(center1,elgvn,eqm,ndipole,ientro,iterld,iprint)
  use chemsol, only : ef_ld,vlgvn_f
  implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000)
      parameter (mxatm=500)
      parameter (mxcenter=50)
      PARAMETER (ONE=1.0d0, ZERO=0.d0)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
           n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm), &
           q_gas(mxatm),q_mp2(mxatm)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn), &
           da(3,mxlgvn)
      common /pcgrid/ drg,rg,dxp0(3), &
           rg_inner,drg_inner
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common/scrat5/  rz1(mxlgvn),rz_vdw(mxlgvn), iz(mxlgvn)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /atom_fs/ atomfs(mxatm) 
      common /outer_surf/ isd (mxlgvn)
      common /lra/ clgvn, slgvn

      !  input vars
      real*8 center1(3)
      ! local vars
      integer*2 isd 
      ! for writing binary data
      integer :: irec

      real*8 elgvn,elgvna
      real*8 fs, gri_sp, efn, fma
      real*8 vdwsur(mxatm)
      character*1 dash(72)
      real(8), dimension(3) :: vlgvn_result
      real(8) :: ea = 0.0d0 ! variable for holding ea from vlgvn
      data dash/72*'-'/

      elgvn=0.d0 !total lgvn energy
!     ephil1 and ephil2 are defined in readopt, sres is surface for 
!     large fields.

      sres= 0.0d0 !?
      !elgvn=0.d0                ! variable set to 0 again
      fsurfa(ientro)=0.d0        ! ientro = 1 at first vlgvn call
      evdwl(ientro)=0.d0         
      ndipole=0
      idum = 1
      tds = 0.0d0               ! initialize tds to 0, should be done in lgvnx
      ! after this is called, ndipole is set
      ! n_inner = 0 before call, 741 after call
      call gen_gridx(center1,ndipole,ientro,0,iprint)

! --   Cartesian coordinates of point dipoles {Angstrom} are stored
!      (it is unclear at present why 2 different variables (xl, xd)
!      were used for coordinates of dipoles in the old pdld code.)
      do i=1,ndipole
      xd(1,i) = xl(1,i)
      xd(2,i) = xl(2,i)
      xd(3,i) = xl(3,i)
      end do


!      Calculate magnitudes of langevin dipoles and noniterative solvation
!      energy from the electric field scaled by the distance-dependent 
!      dielectric constant.

! --   Get electric field at the positions of the dipoles (da(3,mxlgvn))
!      call ef_ld(ndipole,1)

!      inquire (iolength=irec) xw
      ! print the solute cartesian coordinates and charge
      ! open(file="atoms.txt",unit=46,status="replace")
      ! do i=1,n_reg1
      !    write(46,*) xw(1,i),xw(2,i),xw(3,i)
      ! end do
      ! close(46)
      ! open(file="atom-charges.txt",unit=46,status="replace")
      ! do i=1,n_reg1
      !    write(46,*) q(i)
      ! end do
      ! close(46)
      ! ! print the langevin dipoles cartesian coordinates
      ! open(file="dipoles.txt",unit=47,status="replace")
      ! do i=1,ndipole
      !    write(47,*) xl(1,i),xl(2,i),xl(3,i)
      ! end do
      ! close(47)
      da = ef_ld(xw,q,n_reg1,xl,ndipole,1)
      ! open(file="da.txt",unit=48,status="replace")
      ! do i=1,ndipole
      !    write(48,*) da(1,i),da(2,i),da(3,i)
      ! end do
      ! close(48)
      
      efn_max=-10.0d0
!      elgvn = 0.d0
      do i = 1,n_reg1
      vdwsur(i) = 0.d0
      end do
      
      do i=1,ndipole

         gri_sp=drg_inner
         if(i.gt.n_inner) gri_sp=drg 
         if(i.eq.n_inner+1) elgvn_i = elgvn ! save elgvn_i for inner grid ld energies report

         efn=dsqrt(da(1,i)*da(1,i)+da(2,i)*da(2,i)+da(3,i)*da(3,i))
         if (efn.gt.efn_max) efn_max=efn ! set a new efn_max if a greater one is found
         ! call vlgvn(efn,elgvn,xjunk,fma,gri_sp)
         vlgvn_result = vlgvn_F(efn,gri_sp,slgvn)
         fma = vlgvn_result(1)
         tds = vlgvn_result(2)
         elgvn  =  elgvn + vlgvn_result(3)
         xmua(1,i)=fma*da(1,i)/efn
         xmua(2,i)=fma*da(2,i)/efn
         xmua(3,i)=fma*da(3,i)/efn


         ! --  Calculate hydrophobic surface
         if (i.lt.n_inner) then
            if (rz_vdw(i).le.rzcut) then

               ! --    Calculate elstat. potential at the grid point
               epot = 0.d0
               do k=1,n_reg1
                  rx = xl(1,i) - xw(1,k) 
                  ry = xl(2,i) - xw(2,k) 
                  rz = xl(3,i) - xw(3,k) 
                  rqd = sqrt (rx*rx+ry*ry+rz*rz)
                  epot = epot + q(k)/rqd
               end do

               ! --    Hydrophobic energy - negative surfaces
               !      if (epot.lt.0.d0 .and. epot.gt.-ephil2) then
               !       fs = epot/ephil2 + 1.d0
               !       write(6,'(i5,2f10.5," > - ephil2, fs = ",f10.5)') iz(i),
               !    *  rz_vdw(i), epot, fs
               !      end if
               !     
               if (epot.lt.0.d0) epot=-epot
               ! --    Positive surfaces and the surface for vdw term 
               if (epot.le.ephil1) then
                  fs=1.0d0
                  !        write(6,'(i5,2f10.5," < ephil1, fs = ",f10.5)') iz(i),
                  !    *   rz_vdw(i), epot, fs
               end if
               if((epot.gt.ephil1).and.(epot.le.ephil2)) then
                  fs=1.d0-((epot-ephil1)/(ephil2-ephil1))*(1.0-sres)
                  !        write(6,'(i5,2f10.5," < ephil2, fs = ",f10.5)') iz(i),
                  !    *   rz_vdw(i), epot, fs
               elseif(epot.gt.ephil2) then
                  fs=sres
                  !        write(6,'(i5,2f10.5," > ephil2")') i, rz_vdw(i), epot
               endif
               !      endif
               atomfs(iz(i)) = atomfs(iz(i))+fs
               fsurfa(ientro)=fsurfa(ientro)+fs
               vdwsur(iz(i)) = vdwsur(iz(i)) + 1.d0
            endif
         endif
         ! #GOTO 10
      enddo

!     Use atom polarizabilities (vdwc6) to calculate vdW part
!     of the solvation enthalpy.

      evdwl(ientro) = 0.d0
      do k=1,n_reg1
      evdwl(ientro) = evdwl(ientro) + vdwc6(iacw(k))*vdwsur(k)
      end do
      evdwl(ientro) = vdwsl*evdwl(ientro)

      elgvn = clgvn * elgvn
      elgvn_i = clgvn * elgvn_i
!     write(6,'(/,"Maximum field = ",f10.4,/)') efn_max
!------------------------------------------------------------
      write(6,1001) ndipole,elgvn, elgvn_i ! need to get here
!-------------------------------------------------------------
      if (iterld .eq. 0) return

!      Calculation of the initial configuration of Langevin dipoles. 
!      (0-th step of the iterative calculation of dipole-dipole interactions).

      da = ef_ld(xw,q,n_reg1,xl,ndipole,0)

      elgvna=0.0
      do 20 i=1,ndipole
      gri_sp=drg_inner
      if(i.gt.n_inner) gri_sp=drg
      efn=dsqrt(da(1,i)*da(1,i)+da(2,i)*da(2,i)+da(3,i)*da(3,i))
      ! call vlgvn(efn,elgvna,dumm,fma,gri_sp)
      vlgvn_result = vlgvn_F(efn,gri_sp,slgvn)
      fma = vlgvn_result(1)
      tds = vlgvn_result(2)
      elgvn  = elgvn + vlgvn_result(3)
      
! --  Dipole moments of point (langevin) dipoles are oriented along
!     the field from solute for inner grid and outer surface dipoles. 
!     They are oriented randomly for odd nonsurface outer grid 
!     dipoles. (Note that this is implemented 
!     for the zero iteration step only. For iterative lgvn,
!     lgvn dipoles are calculated in subroutine mu_mu_l. Only those
!     of lgvn dipoles that lie along outer surface are constrained 
!     to be proportional to the solute field.) 

      if(i.le.n_inner .or. isd (i).eq.1 .or. mod(i,3).eq.0) then 
      ddd=3.0d0/(sqrt(rz1(i))+2.0d0)
      xmua(1,i)=ddd*fma*da(1,i)/efn
      xmua(2,i)=ddd*fma*da(2,i)/efn
      xmua(3,i)=ddd*fma*da(3,i)/efn
      go to 20

!  -- The remaining outer grid dipoles are placed randomly.
      else
      dddx=ran2(idum)
      dddy=ran2(idum)
      dddz=ran2(idum)
      ddd=sqrt(dddx*dddx+dddy*dddy+dddz*dddz)
      xmua(1,i)=dddx*fma/ddd
      xmua(2,i)=dddy*fma/ddd
      xmua(3,i)=dddz*fma/ddd
      end if
20    continue

      return
!......................................................................
1001  format(' Total number of Langevin dipoles      : ',i10/ &
           ' Noniterative lgvn energy (total,inner): ',2f10.3/)
      end
