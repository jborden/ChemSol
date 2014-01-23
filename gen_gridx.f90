      subroutine gen_gridx (center1,nld,ientro,iflag,iprint)
      implicit Real*8 (a-h,o-z)
!     generate grid point for langevin dipole
!     center1...center of the grid (input)
!     nld...number of grid points (output)
!

      parameter (mxatm=500)
      parameter (mxlgvn=10000)
      parameter (mxcenter=50)

      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
           n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm) &
           ,q_gas(mxatm),q_mp2(mxatm)
      common /born/ rg_reg1, rgim
      common /pcgrid/ drg,rg,dxp0(3), &
           rg_inner,drg_inner
!:::  center
      common /pcgmcntr/ pcenter(3)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn), &
           da(3,mxlgvn)
      common/scrat5/  rz1(mxlgvn),rz_vdw(mxlgvn), iz(mxlgvn)
!::   for surface and vdw calculation
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /aname/ atom(mxatm),molname,ssname
      common /outer_surf/ isd (mxlgvn)
      common /volume/ nvol(mxcenter)
      character*8 atom
      character*13 molname
      character*4 ssname

!...  input vars
      integer i0,i1,nld,iflag
      integer*2 isd
      real*8 center1(3)
 
!...  local vars
      integer i,index,ii,jj,kk,no,n_out1
      integer ipro,jpro,kpro,igrid,jgrid,kgrid
      integer icube(0:133,0:133,0:133)
!     integer limit_inner,limit_outer,mid_inner,mid_outer,idone
      real*8 rg2,rgi2,rgim2,drg2,drgi2,drg_i,drg_inner_i
      real*8 rp2(mxatm)
!     save limit_inner, limit_outer,mid_inner,mid_outer, 
!    $     rg2,rgi2,rgim2,drg2,drgi2,drg_i,drg_inner_i,idone
      real*8 ri,rj,rk,xloca,yloca,zloca,vdw_in,vdw_out
      real*8 xp(3),d2,xdrg,r6,d6,r3,d3
!     data idone/0/
!.......................................................................
!     Set up parameters. 
!     rg is a grid radius (outer), rgi is the same for the inner grid.
!     drg is outer grid spacing (3A), rgim is a trim parameter
!     for the inner grid.
      nld=0
      ns=0
      i0=1
      i1=n_reg1
      rg2=rg*rg
      drg2=drg*drg
      drgi2=drg_inner*drg_inner
      rgi2=rg_inner*rg_inner
      rgim2=rgim*rgim
      drg_inner_i=1.d0/drg_inner
      drg_i=1.0d0/drg
      do i=1,n_reg1
      rp2(i)=rpi(i)*rpi(i)
      end do

      limit_inner=int(rg_inner*2/drg_inner)
      limit_outer=int(rg*2/drg)
!*
!*     [make both limits into an odd number]
!*     
      if(mod(limit_inner,2).eq.0) limit_inner=limit_inner+1
      if(mod(limit_outer,2).eq.0) limit_outer=limit_outer+1
      mid_inner=int(limit_inner/2.d0+0.5d0)
      mid_outer=int(limit_outer/2.d0+0.5d0)

!......................................................
!     build inner grid, if used.
!......................................................

      if(rg_inner.gt.0.d0) then ! build inner grid
         do  ii=0,limit_inner
            do jj=0,limit_inner
               do  kk=0,limit_inner
                  icube(ii,jj,kk)=0
                  ri=ii-mid_inner
                  rj=jj-mid_inner
                  rk=kk-mid_inner
                  d2=(ri*ri+rj*rj+rk*rk)*drgi2 !radius from center
                  if(d2.le.rgi2)icube(ii,jj,kk)=1   !trim to a sphere
               enddo
            enddo
         end do

!     Make cavity 
!     Put solute on the grid nad find the nearest grid points
!     Out of these, grid points closer than rp will be removed.
!     The number of points removed (nvol) is proportional to the solute
!     volume.

      nvol(ientro) = 0

      do i=i0,i1 !loop over solute atoms
         ipro=int(mid_inner+(xw(1,i)-center1(1))*drg_inner_i)
         jpro=int(mid_inner+(xw(2,i)-center1(2))*drg_inner_i)
         kpro=int(mid_inner+(xw(3,i)-center1(3))*drg_inner_i)
         if(ipro.ge.0.and.ipro.le.limit_inner.and. &
              jpro.ge.0.and.jpro.le.limit_inner.and. &
              kpro.ge.0.and.kpro.le.limit_inner) then !within grid/dimensi
!           [check distances to 27 closest points around the
!            nearest grid point including itself]
            do ii=1,9
               igrid=ipro+ii-5
               do jj=1,9
                  jgrid=jpro+jj-5
                  do kk=1,9
                     kgrid=kpro+kk-5
                     if(icube(igrid,jgrid,kgrid).eq.1)then
                        xloca=(igrid-mid_inner)*drg_inner+center1(1)
                        yloca=(jgrid-mid_inner)*drg_inner+center1(2)
                        zloca=(kgrid-mid_inner)*drg_inner+center1(3)
                        ri=xloca-xw(1,i)
                        rj=yloca-xw(2,i)
                        rk=zloca-xw(3,i)
                        d2=ri*ri+rj*rj+rk*rk !distance from grid point
                        if(d2.lt.rp2(i)) then
                           icube(igrid,jgrid,kgrid)=0
                           nvol(ientro) = nvol(ientro) + 1
                        end if
                     endif
                  end do
               end do
            end do
         endif
      end do
 
!     Reshape grid sphere (outer region) to a solute envelope.
      do ii=0,limit_inner
         do jj=0,limit_inner
            do kk=0,limit_inner
               if(icube(ii,jj,kk).ne.0) then
                  xp(1)=center1(1)+(ii-mid_inner)*drg_inner
                  xp(2)=center1(2)+(jj-mid_inner)*drg_inner
                  xp(3)=center1(3)+(kk-mid_inner)*drg_inner
                  d2_min=10000.d0
                  !                 Get distance to a closest VdW boundary.
                  do i=i0,i1
                     ri=xp(1)-xw(1,i)
                     rj=xp(2)-xw(2,i)
                     rk=xp(3)-xw(3,i)
                     d2=ri*ri+rj*rj+rk*rk
                     d2=sqrt(d2)
                     d2=d2-rpi(i)
                     if (d2.lt.0.d0) then
                        write (6,'("Wrongly placed grid point in gen_gridx!")')
                        stop
                     elseif(d2.le.d2_min) then
                        d2_min=d2
                     end if
                  end do
                  if(d2_min.le.rgim) then !less than rgim
                     nld=nld+1
                     !......................................................................2
                     if(nld.gt.mxlgvn)then
                        nld=nld-1
                        write(6,'("Exceeds current Dipole grid limits use smaller rg in option file!")')
                        stop
                     endif
                     isd (nld)=0
                     xl(1,nld)=xp(1)
                     xl(2,nld)=xp(2)
                     xl(3,nld)=xp(3)
                  endif
               endif
            end do
         end do
      end do
   endif

      n_inner=nld
      write(6,1003) nld 
1003  format(/' No of inner grid dipoles              : ',i10)
!......................................................
!     build outer grid (3.0 spacing)
!......................................................
 
      do ii=0,limit_outer
         do jj=0,limit_outer
            do kk=0,limit_outer
               icube(ii,jj,kk)=0
               ri=ii-mid_outer
               rj=jj-mid_outer
               rk=kk-mid_outer
               d2=(ri*ri+rj*rj+rk*rk)*drg2
               if(d2.le.rg2) icube(ii,jj,kk)=1 !within rg
            end do
         end do
      end do
 
!     Make cavity 
!     Put solute on the grid and find the nearest grid points.
!     Out of these, grid points closer than rp will be removed.
!     Also, grid points coordinates of which are less than 2A from
!     the inner grid  points are removed.

      do i=i0,i1 
!        [find distance to nearest solute atom]
         ipro=int(mid_outer+(xw(1,i)-center1(1))*drg_i)
         jpro=int(mid_outer+(xw(2,i)-center1(2))*drg_i)
         kpro=int(mid_outer+(xw(3,i)-center1(3))*drg_i)
         if(ipro.ge.0.and.ipro.le.limit_outer.and. &
              jpro.ge.0.and.jpro.le.limit_outer.and. &
              kpro.ge.0.and.kpro.le.limit_outer) then !within grid/dimensi

            do ii=0,4
               igrid=ipro+ii-2
               do jj=0,4
                  jgrid=jpro+jj-2
                  do kk=0,4
                     kgrid=kpro+kk-2
                     if(icube(igrid,jgrid,kgrid).eq.1)then
                        xloca=(igrid-mid_outer)*drg+center1(1)
                        yloca=(jgrid-mid_outer)*drg+center1(2)
                        zloca=(kgrid-mid_outer)*drg+center1(3)
                        ri=xloca-xw(1,i)
                        rj=yloca-xw(2,i)
                        rk=zloca-xw(3,i)
                        d2=ri*ri+rj*rj+rk*rk !distance from grid point
                        if(d2.lt.rp2(i)) icube(igrid,jgrid,kgrid)=0
                     endif
                  end do
               end do
            end do
         end if
      end do

      no=0
      ns=0
      efmin1=0.0015d0
      efmin2=0.0020d0
!     efmin1=0.0020d0
!     efmin2=0.0024d0

      do 41 ii=0,limit_outer
      do 41 jj=0,limit_outer
      do 41 kk=0,limit_outer
         if(icube(ii,jj,kk).ne.0) then
            no=no+1
            xp(1)=center1(1)+(ii-mid_outer)*drg
            xp(2)=center1(2)+(jj-mid_outer)*drg
            xp(3)=center1(3)+(kk-mid_outer)*drg
!           [loop over inner grid points]
!           adjust spacing between the inner and outer grid
            do index=1,n_inner
               ri=xl(1,index)-xp(1)
               rj=xl(2,index)-xp(2)
               rk=xl(3,index)-xp(3)
               d2=ri*ri+rj*rj+rk*rk
               d2=sqrt(d2)
               if(d2.lt.((drg_inner+drg)/2.-0.001)) goto 41
!              if(d2.lt.drg-0.001) goto 41
            end do

! -- Remove grid points in the areas with small electric field from
!    the solute. Use distance-dependent dielectric.
! -- For +1 charged ion, the ef>0.002  e/A**2 threshold results in the
!    grid radius of 14.4 A.
! -- Finaly, gridpoints with 0.002>ef>0.0024 will be marked in the
!    list isd (mxlgvn). Lgvn dipoles at these grid points will be kept
!    constant in future calculations (surface constrained lgvn dipoles).

            efx=0.d0
            efy=0.d0
            efz=0.d0

            do i=1,n_reg1
               ri=xp(1)-xw(1,i)      
               rj=xp(2)-xw(2,i)      
               rk=xp(3)-xw(3,i)      
               d2=ri*ri+rj*rj+rk*rk
               d1=sqrt(d2)
               ddd=1.7d0/sqrt(d1+2.d0)
               d3=d1 * d2
               qr=ddd*q(i)/d3
               efx=efx + qr*ri
               efy=efy + qr*rj
               efz=efz + qr*rk
            end do
            efnorm=sqrt(efx**2+efy**2+efz**2)
            if (efnorm.lt.efmin1) goto 41
            nld=nld+1
            if(nld.gt.mxlgvn)then
               nld=nld-1
               write(6,'("Exceeds current Dipole grid limits use smaller rg in option file!")')
               stop
            endif

            if(efnorm.lt.efmin2) then
            isd (nld)=1
            ns = ns + 1
            else
            isd (nld)=0
            end if

            xl(1,nld)=xp(1)
            xl(2,nld)=xp(2)
            xl(3,nld)=xp(3)
         endif
 41   continue

      n_out1=nld-n_inner
      write(6,1001) ns
1001  format(' No of surface constrained dipoles     : ',i10) 

!     Change the positions of outer-grid dipoles randomly in order to 
!     eliminate artifacts related to the regularity of the cubic 
!     grid.
!
!     do i = n_inner+1, nld
!       xl(1,i) = xl(1,i) + ran2(idum)/2.d0 -0.5d0
!       xl(2,i) = xl(2,i) + ran2(idum)/2.d0 -0.5d0
!       xl(3,i) = xl(3,i) + ran2(idum)/2.d0 -0.5d0
!     end do

!     Collect distances of grid dipoles to solute nuclei (rz1) and 
!     solute boundary (rz_vdw). Note that rz1 is a SQUARE of the dist.
      do i=1,nld
         rz1(i)=10000.  
         d2_min=10000.
         do index=1,n_reg1 
            ri=xl(1,i)-xw(1,index)
            rj=xl(2,i)-xw(2,index)
            rk=xl(3,i)-xw(3,index)
            d2=ri*ri+rj*rj+rk*rk
            rz1(i)=dmin1(rz1(i),d2)
            d2=sqrt(d2)
            d2=d2-rpi(index)
            if (d2.lt.0.d0) then
            write (6,'("Grid point, n = ",i6," within Vdw radius of atom: ",i5,f10.3," is less than rp!")') i,index,d2    
            stop
            elseif(d2.le.d2_min) then
            d2_min=d2
            iz(i)=index
            rz_vdw(i)=d2
            end if
         enddo
!        write(6,'(i7,i6,f10.3)') i, iz(i), rz_vdw(i) 
      enddo

! --  VdW energy (9-6 formula)
!     The minimum of the VdW curve is at the rp distance.
!     London coeficients (vdwC6) are atom dependent, but they are not 
!      hybridization-dependent. (e.g, (vdwc6(C(sp3)) = vdwc6(C(sp2)). 

      evdwl(ientro)=0.0
      vdw_in=0.d0
      vdw_out=0.d0

!###################SKIPPED-BEGIN##########################
! -- Calculate vdW energy elsewhere
 
      if (.false.) then

      if (iprint.eq.1) then
         write(6,'("drg_i,drg,vdwsl= ",3f8.2)') drg_inner, drg, vdwsl
         write(6,'("n_inner, nld =  ",2i8)') n_inner, nld
         write(6,'(/,"Contribution to VdW from individual atoms",/)')
      end if

      do index=i0,i1
         subvdw = vdw_in + vdw_out
         r3 = rp2(index)*rpi(index)
         r6 = r3 * r3
         c6 = vdwc6(iacw(index))/rp2(index)
         xdrg=(drg_inner/3.d0)**3
         do i=1,n_inner
            ri=xl(1,i)-xw(1,index)
            rj=xl(2,i)-xw(2,index)
            rk=xl(3,i)-xw(3,index)
            d2=ri*ri+rj*rj+rk*rk
            d3=d2*sqrt(d2)
            d6=d2*d2*d2
            vdw_in=vdw_in+vdwsl*c6*(2.0*r6*r3/d6/d3-3.0*r6/d6)*xdrg
         end do

         xdrg=(drg/3.d0)**3
         do i=n_inner+1,nld
            ri=xl(1,i)-xw(1,index)
            rj=xl(2,i)-xw(2,index)
            rk=xl(3,i)-xw(3,index)
            d2=ri*ri+rj*rj+rk*rk
            d3=d2*sqrt(d2)
            d6=d2*d2*d2
            vdw_out=vdw_out+vdwsl*c6*(2.0*r6*r3/d6/d3-3.0*r6/d6)*xdrg
         end do
         subvdw=vdw_in+vdw_out-subvdw
         if(iprint.eq.1) then
            write(6,'(1x,i4,2x,a8,f10.5)')index,atom(index),subvdw
         end if
      end do

      evdwl(ientro)=vdw_in+vdw_out
      write(6,'(/,"Total vdw: in, out, tot: ", 3f10.5)') vdw_in, &
           vdw_out, evdwl(ientro)

   end if
     
   return
end subroutine gen_gridx
