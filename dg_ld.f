c--------------------------------------------------------------------
      subroutine dg_ld (iterld,iprint)
      implicit Real*8 (a-h,o-z)

      parameter (mxcenter=50)
      PARAMETER (MXATM=500)
      common /aname/ atom(mxatm),molname,ssname
      common /pctimes/ ndxp,itl,itp
      common /pcgmcntr/ pcenter(3)
      common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82),
     $      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)
     :             ,q_gas(mxatm),q_mp2(mxatm)
      common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /atom_fs/ atomfs(mxatm)
      character*8 atom
      character*13 molname
      character*4 ssname

c:::  local vars
      integer i,j
      real*8 center_new(3)
      character*1 dash(72)
      data dash/72*'-'/
c.......................................................................
      esum = 0.d0
      evqdq = 0.d0
      ecor=0.d0
      do i=1,n_reg1
      atomfs(i)=0.0d0
      end do

c      elgvn.......noniterative lgvn energy (using distance-dependent dielectric)
c      elgvni.....lgvn energy from the 0-th iteration
c      elgwa......iterative lgvn energy 
c      evqdq......solute relaxation energy calculated from PCM charges
c      ecor.......correction for electron correlation effects         

      do 10 i=1,ndxp
         call ran_shift(i,pcenter,center_new)
c        Set grid origin for the printout of dipoles into an xmol input file.
c        center_new(1)=0.0
c        center_new(2)=0.0
c        center_new(3)=0.0
         call lgvnx(center_new,elgvn,eqm,ndipole,i,iterld,iprint)
         esum = esum + elgvn

         if (iterld.eq.1) then
         call sci_lgvn(elgwa,elgvni,tds,ndipole,i)
         temp_elgvn(i)=elgwa
         tdsl(i) = tds
         else
         elgwa= 0.d0
         temp_elgvn(i)=elgvn
         tdsl(i) = 0.d0
         end if

         tdsw_a = phobsl * fsurfa(i)
         if (iprint.eq.1) then
         write(6,'(/," Averaged atomic contributions to E_phobic:",/)')
         do j=1,n_reg1
         write(6,112) j, atom(j), phobsl*atomfs(j)/float(i)
         end do
         end if
         if (iterld.eq.1) write(6,105) elgwa,evdwl(i),tdsw_a,-tds+tdsw_a
         if (iterld.eq.0) write(6,105) elgvn,evdwl(i),tdsw_a,-tds+tdsw_a
         write(6,102) dash

         call vatom (cor,vqdq,ndipole,iterld,iprint)
         evqdq = evqdq + vqdq
         ecor = ecor + cor
 10   continue

      elgvn = esum/dble(ndxp)
      evqdq = - evqdq/dble(ndxp)
      ecor  = ecor/dble(ndxp)
      call elgvn_ave(iterld,ndxp)
      call vbornx(ndipole,ebw,center_new)

      return
c.......................................................................
 102  format(1x,72a1)
 105  format(/1x,'Elgvn = ',f8.3,'   Evdw = ',f8.3,
     $       '   -TdS_phobic = ',f8.3,'   -TdS_total = ',f8.3,/)
 112  format(i3,5x,a8,f9.3) 

      end
