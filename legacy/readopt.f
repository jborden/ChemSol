
c--------------------------------------------------------------------
      subroutine readopt (iterld)
      implicit Real*8 (a-h, o-z)

      PARAMETER (MXATM=500)
      parameter (mxcenter=50)
      character*1 dash(72)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82),
     :            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)
     :             ,q_gas(mxatm),q_mp2(mxatm)
      common /born/ rg_reg1, rgim
      common /aname/ atom(mxatm),molname,ssname
      common /pcgmcntr/ pcenter(3)
      common /pcgrid/ drg,rg,dxp0(3),
     :              rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /pctimes/ ndxp,itl,itp
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /lra/ clgvn, slgvn
      character*8 atom
      character*13 molname
      character*4 ssname
      character*2 rpinp
      data dash/72*'-'/
c....................................................................
C     For the adjustment of VdW radii (rp) and London coef (vdwc6) 
C     additional input from vdw.par file is added here:
      do i=1,82
      vdwc6(i) = 0.
      end do
      read(44,1202) iterld, ndxp, dxp0(1), clgvn, slgvn,tds0
      read(44,1200) srp, rp(6), rp(7), rp(9), rp(10)
      read(44,1200) rp(13), rp(14), rp(15), rp(22), rp(24)
      read(44,1200) rp(16), rp(25)
      read(44,1200) vdwc6(1), vdwc6(6), vdwc6(9), vdwc6(13),vdwc6(22)
      read(44,1200) vdwc6(24), vdwc6(7)
      read(44,1200) vdwsl, phobsl, ephil1,ephil2, rzcut
      phobsl=phobsl/10.d0
      vdwsl=-vdwsl/10.d0

c      Define the remaining London coef 
      do i=1, n_reg1
      if (vdwc6(iacw(i)) .eq.0.d0 ) then
        if(vdwc6(iacw(i)-1).ne.0) then
        vdwc6(iacw(i)) = vdwc6(iacw(i)-1)
        else
        vdwc6(iacw(i)) = vdwc6(iacw(i)-2)
        end if
      end if
      if (vdwc6(iacw(i)) .eq.0.d0 ) then
      if (iacw(i).le.17) vdwc6(iacw(i))=vdwc6(6)
      if (iacw(i).gt.17) vdwc6(iacw(i))=vdwc6(22)
      end if
      end do

c     Calculate rp(H) as a linear function of the rp of the atom to        
c     which hydrogen is covalently bonded. Note that the coefficient of 
c     this function (srp) differs for 1st and 2nd row atoms.
c     In addition, for inorganic oxygen (iacw=15)
c     a separate rp(H) is used (it is read from vdw.par file).
      do 60 i=1,n_reg1
      if(iacw(i).ne.1 .and.iacw(i).ne.2) then
      rpi(i) = rp(iacw(i))
      else
      r2min = 100.
      jmin = 21
      do 50 j=1,n_reg1
        if(j.eq.i) go to 50
        r2=(xw(1,i)-xw(1,j))**2+(xw(2,i)-xw(2,j))**2
     $      +(xw(3,i)-xw(3,j))**2
        if(r2.lt.r2min) then
          r2min=r2
          jmin=j
        end if
50    continue
         if((iacw(jmin)).eq.15) then
         iacw(i) = 2
         rpi(i) = rp(iacw(i))
         vdwc6(iacw(i)) = vdwc6(iacw(i)-1)
         else
         if(iacw(jmin).lt.18) rpi(i) = srp * rp(iacw(jmin))
         if(iacw(jmin).ge.18) rpi(i) = (srp-0.1)*rp(iacw(jmin))
         end if
      end if
60    continue

C --  Finally, allow for the change of rp parameter of any atom
C     without changing its atom type. A new rp is specified in
C     the end of the input file.
      read (45,'(a2)') rpinp
      if (rpinp.eq.'rp') then
      read (45,*) nrp
      do j=1,nrp
      read (45,*) i, rpi(i)
      end do
      end if

      if (iprint.eq.1) then
      write(6,'(/,"LD parameters: ",/)')
      write(6,1202) iterld, ndxp, dxp0(1), clgvn, slgvn
      write(6,1200) vdwsl, phobsl, ephil1, ephil2, rzcut
      end if

      write(6,95)
      do i=1,n_reg1
      write (6,102) latom(i),(xw(j,i),j=1,3),q(i),q_gas(i),
     $              iacw(i), rpi(i),vdwc6(iacw(i))
      enddo
      write(6,103) dash

 95   format(//,1x,'atom #',5x,'x',8x
     $,'y',8x,'z',3x,'Q_pcm',2x,'Q_gas',1x,'atom_type',
     $3x,'rp',4x,'VdWC6'/)
 102  format (i6,3f9.3,2f6.2,i6,7x,f4.2,4x,f4.2)
 103  format (/1x,72a1/)
1200  format (5f6.3)
1202  format (i2,i3,4f6.3)

c --  Write input file for Polaris (This file is used by x_prep
c     to create entry into amino-acid library and to create pdb file.
c      call getenv('XPOLOUT',pname)
c      open (36, file=pname)
c      write(36,'(a13)') molname
c      do i=1,n_reg1
c      nbond=0
c      do jj=1,10
c      nb(jj)=0.
c      end do
c      do j=1,n_reg1
c        d1=(xw(1,i)-xw(1,j))**2+(xw(2,i)-xw(2,j))**2
c     $      +(xw(3,i)-xw(3,j))**2
c        d1=sqrt(d1)
c        if (i.ne.j.and.(d1.lt.1.8).and.(iacw(i).ne.1.or.
c     $      iacw(j).ne.1)) then
c        nbond=nbond+1
c        nb(nbond)=j
c        end if
c      end do
c      write(36,1211) i, atom(i), molname, xw(1,i), xw(2,i),
c     $               xw(3,i), q(i), nbond, (nb(k),k=1,nbond)
c      end do
c1211  format(i3,1x,A4,2x,a3,4f9.3, 10i3)
c      close (36)

      pcenter(1)=0.
      pcenter(2)=0.
      pcenter(3)=0.

      do i=1,n_reg1
         pcenter(1)=pcenter(1)+xw(1,i)/n_reg1
         pcenter(2)=pcenter(2)+xw(2,i)/n_reg1
         pcenter(3)=pcenter(3)+xw(3,i)/n_reg1
      enddo

C -- Find maximal radius of the solute wrt grid center.
      dmax=0.d0
      do i=1,n_reg1
      rg_reg1=(pcenter(1)-xw(1,i))**2 + (pcenter(2)-xw(2,i))**2 +
     $        (pcenter(3)-xw(3,i))**2
      if (rg_reg1.gt.dmax) dmax=rg_reg1
      end do
      rg_reg1=sqrt(dmax)+2.4d0
C -- New grid radii are measured wrt solute surface. This ensures
C    that enough space is attributed to the grid for large molecules.
C    For actual choice of grid extension, other criteria are applied:
C    The distance from the VdW surface (inner, 1A grid) and magnitude of 
C    the field at the grid point (outer, 3A grid) - see gen_gridx subroutine.
      rg=rg_reg1+rg
      rg_inner=rg_reg1+rgim+2.d0

      if (iprint.eq.1) then
      write(6,106) pcenter
      write(6,107) rg_reg1
      end if

      if (drg.lt.3.) then
         write(6,996) 
         stop
      elseif (ndxp.gt.mxcenter) then
         write(6,997) mxcenter
         stop
      endif
      
 106  format(/,1x,'Center of the cubic grid is defined as centroid of',
     s     ' solute atoms',//,20x,'center -',3f8.3)
 107  format(/,'Radius of the centroid: ',f6.2,//)
 996  FORMAT(//' PROGRAM EXIT: DRG SHOULD BE AT LEAST 3.O A'//)
 997  FORMAT(//' PROGRAM EXIT: MAXIMUN # OF GRIDS IS',I3//)
      
      return
      end
