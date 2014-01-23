      subroutine pairlistw(nd,dipx)
!     subroutine pairlistw(nd,dipx,jp3)
! --  The grid-point pairs are put here on the three separate lists
!     according the distance between them. These lists are subsequently
!     used in the calculation of dipole-dipole interactions.
!     Note that interactions between nearest neighbours (1A grid)
!     and long-range interactions (>20A) are excluded.

      implicit Real*8 (a-h,o-z)
      parameter (mxatm=500)
      parameter (mxpair=5000000)
      parameter (mxpair2=10000000)
!     parameter (mxpair3=1800000)
      parameter (mxlgvn=10000)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
           n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm), &
           q_gas(mxatm),q_mp2(mxatm)
      common /pcgrid/ drg,rg,dxp0(3), &
           rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn), &
           jp2(mxpair2),ip3(0:mxlgvn)
      common /outer_surf/ isd (mxlgvn)

      integer*2 jp,jp2,isd
      integer   ip,ip2,ip3
      integer nd
!     integer*2 jp3(mxpair3)
      real*8 dipx(3,mxlgvn)

!:::  local vars
      integer i,j
      integer npair,npair2,npair3,ntotl
      real*8 d2,rdcut2,r1,r2,r3,out_cut2
      integer mask1
      real*8 rcoff,rcoff2,out_cut
!....................................................................
      rdcut2=rdcutl*rdcutl
      out_cut2=out_cut*out_cut

      npair=0
      ip(0)=0

      npair2=0
      ip2(0)=0

      npair3=0
      ip3(0)=0

      ntotl=0

      do 10 i=1,nd-1
         do 15 j=i+1,nd
            if(isd(i).eq.1.and.isd(j).eq.1) goto 15
            r1=dipx(1,i)-dipx(1,j)
            r2=dipx(2,i)-dipx(2,j)
            r3=dipx(3,i)-dipx(3,j)
            d2=r1*r1+r2*r2+r3*r3
            rcoff=2.5
!  -- Switch on/off nearest neighbors for outer grid:
!           if(i.gt.n_inner.and.j.gt.n_inner) rcoff=drg+0.2
!  -- Switch on/off interaction of nearest neighbors belonging to 
!     the inner and outer grid:
!           if(i.le.n_inner.and.j.gt.n_inner) rcoff=1.9
            rcoff2=rcoff*rcoff
            mask1=d2/rcoff2
            if(mask1.ge.1) then
               ntotl=ntotl+1
               if(d2.lt.rdcut2) then
                  npair=npair+1
                  if (npair.gt.mxpair) then
                     write(6,1000)npair,mxpair
                     stop 
                  endif
                  jp(npair)=j
               else if(d2.ge.rdcut2.and.d2.lt.out_cut2)then
                  npair2=npair2+1
                  if (npair2.gt.mxpair2) then
                     write(6,1002) npair2,mxpair2
                     stop 
                  endif
                  jp2(npair2)=j
!              else if(d2.ge.out_cut2)then
!                 npair3=npair3+1
!                 if (npair3.gt.mxpair3) then
!                    write(6,1003)npair3,mxpair3
!                    stop 
!                 endif
!                 jp3(npair3)=j
               endif
            endif
! #GOTO 15            
 15      continue
         ip(i)=npair
         ip2(i)=npair2
         ip3(i)=npair3
! #GOTO 10
 10   continue

      write(6,1001) ntotl,npair,rdcutl,npair2,out_cut, &
           npair3,out_cut

      return
!....................................................................
 1000 format(/,1x,i8,' exceeds program mxpair dimension - ',i8)
 1002 format(/,1x,i8,' exceeds program mxpair2 dimension - ',i8)
!1003 format(/,1x,i8,' exceeds program mxpair3 dimension - ',i8,/
!   $         1x,'please use a smaller RDCUTL'//)
 1001 format(/1x,'dipole pair list'/ &
           /,1x,'       total possible pairs  - ',i10,' out of ', &
           '   2.500 A cutoff', &
           /,1x,'       pairs selected        - ',i10,' within ', &
           f8.3,' A cutoff', &
           /,1x,'       pairs selected (long) - ',i10,' within ', &
           f8.3,' A cutoff', &
           /,1x,'       pairs selected (long) - ',i10,' out of ', &
           f8.3,' A cutoff')
   end subroutine pairlistw
