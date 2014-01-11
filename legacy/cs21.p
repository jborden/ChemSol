      implicit Real*8 (a-h,o-z)
      PARAMETER (MXATM=500 )
      PARAMETER (ONE=1.0d0 , ZERO=0.0d0 )
      parameter (mxcenter=50 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),             n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /aname/ atom(mxatm),molname,ssname
      common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
      common /born/ rg_reg1, rgim
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /pctimes/ ndxp,itl,itp
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /lra/ clgvn, slgvn
      character*8 atom,dumm1
      character*13 molname
      character*4 ssname
      character*256 fname
      character*256 input_file_name
      character*256 vdw_name
      character*3 relax
      logical do_gas
      integer iac_conv (89 ), iacp(mxatm), mass(88 )
      data iac_conv/1,2,              3,4,5,6,9,13 ,16 ,17 ,              18 ,19 ,20 ,21 ,22 ,24 ,25 ,26 ,     27 ,28 ,29 ,30 ,31 ,32 ,33 ,34 ,35 ,36 ,37 ,38 ,39 ,40 ,41 ,42 ,43 ,44 ,     45 ,46 ,47 ,48 ,49 ,50 ,51 ,52 ,53 ,54 ,55 ,56 ,57 ,58 ,59 ,60 ,61 ,62 ,     63 ,64 ,65 ,     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,              66 ,67 ,68 ,69 ,70 ,71 ,72 ,73 ,74 ,75 ,76 ,77 ,78 ,79 ,80 ,     81 ,82 /
      data mass /1, 4,           7,9,11 ,12 ,14 ,16 ,19 ,20 ,           23 ,24 ,27 ,28 ,31 ,32 ,35 ,40 ,           39 ,40 ,45 ,48 ,51 ,52 ,55 ,56 ,59 ,59 ,64 ,65 ,70 ,73 ,75 ,79 ,80 ,84 , 85 ,88 ,89 ,91 ,93 ,96 ,97 ,101 ,103 ,106 ,108 ,112 ,115 ,119 ,122 ,128 ,127 ,131 , 132 ,137 ,139 , 140 ,141 ,144 ,145 ,150 ,152 ,157 ,159 ,163 ,165 ,167 ,169 ,173 ,175 , 178 ,181 ,184 ,186 ,190 ,192 ,195 ,197 ,201 ,204 ,207 ,208 ,209 ,210 ,222 , 223 ,226 /
      call getarg(1,vdw_name)
      open (44 ,file=vdw_name)
      call getarg(2,input_file_name)
      open (45 ,file=input_file_name)
      read (45 ,"(a13)") molname
      read (45 ,*) n_reg1, ngeom
      imp2 = 0
      read(45 ,"(a4)") ssname
      do i=1,n_reg1
      read (45 ,1000 ) atom(i),zan(i),q(i),xw(1,i), xw(2,i), xw(3,i)
      iacp(i)=int(zan(i))
      iacw(i)=iac_conv(iacp(i))
      latom(i)=i
      end do
      amas = 0.0d0 
      do i=1,n_reg1
      amas =amas + dble(mass(int(zan(i))))
      end do
      do_gas =  %FALSE% 
      read (45 ,"(a3)") relax
      if (relax EQUAL "pcm") then
        do_gas =  %TRUE% 
        do i=1,n_reg1
          q_gas(i) = q(i)
        end do
        do i=1,n_reg1
          read (45 ,"(17x,f10.4)") q(i)
        end do
      else
        backspace(45 )
        do i=1,n_reg1
          q_gas(i) = q(i)
        end do
      end if
      if(imp2 EQUAL 1) then
      read (45 ,"(a3)") corr
      do i=1,n_reg1
      read (45 ,"(17x,f10.4)") q_mp2(i)
      end do
      else
      do i=1,n_reg1
      q_mp2(i) = q(i)
      end do
      end if
      qsum1 = 0.0 
      qsum2 = 0.0 
      qsum3 = 0.0 
      do i=1,n_reg1
      qsum1 = qsum1 + q(i)
      qsum2 = qsum2 + q_gas(i)
      qsum3 = qsum3 + q_mp2(i)
      end do
      weight = amas
      TdS_gas = entropy (weight)
 55   format(6x,"***********************************************************"/6x,"                                                           "/6x,"                       CHEMSOL 2.1                         "/6x,"                                                           "/6x,"              Jan Florian and Arieh Warshel                "/6x,"            University of Southern California              "/6x,"                    Los Angeles, 1999                      "/6x,"                                                           "/6x,"***********************************************************"//      )
      write (6,55 )
      write (6,"(/\" Solute                           : \",           5x,a13)")molname
      write (6,"(\" Total charge for gas-phase solute: \",f9.2)") qsum1
      write (6,"(\" Total charge for PCM solute      : \",f9.2)") qsum2
      write (6,"(\" Molecular weight (a.u.)          : \",f9.2)") amas
      write (6,"(\" Gas-phase transl. entropy (298dS): \",f9.2,            \" kcal/mol\")") TdS_gas
      k1 = 0
      k2 = 0
      k3 = 0
      do i=1,10 
      if (abs( abs (qsum1)-float(i)+1.0 ) <  0.012 ) k1=1
      if (abs (abs (qsum2)-float(i)+1.0 ) <  0.012 ) k2=1
      if (abs (abs (qsum3)-float(i)+1.0 ) <  0.012 ) k3=1
      end do
      if (k1 EQUAL 0) stop "Inconsistent gas-phase charges"
      if (k2 EQUAL 0) stop "Inconsistent pcm charges"
      if (k3 EQUAL 0) stop "Inconsistent correlated charges"
      do i=1,n_reg1
      if (iacw(i) >< 1) then
        if (atom(i)(2:2)  EQUAL  "1") iacw(i) = iacw(i)+1
        if (atom(i)(2:2)  EQUAL  "2") iacw(i) = iacw(i)+2
        if (atom(i)(2:2)  EQUAL  "3") iacw(i) = iacw(i)+3
      end if
      if (iacw(i)  EQUAL 0) then
      write(6,"(\"Unknown atom type on input:\",i3,a8,f4.1)") i,atom(i),za
      stop
      end if
      end do
      do 30 j=1,ngeom
      if (j >= 2) then
         write(6,2000 ) j
         read(45 ,"(a4)") ssname
         do i=1,n_reg1
         read (45 ,1000 ) dumm1,dumm2,q(i),xw(1,i), xw(2,i), xw(3,i)
         end do
         read (45 ,"(a3)") relax
         if (relax EQUAL "pcm") then
            do i=1,n_reg1
            read (45 ,"(17x,f10.4)") q_gas(i)
            end do
            else
            do i=1,n_reg1
            q_gas(i) = q(i)
            end do
         end if
         if (imp2 EQUAL 1) then
            read (45 ,"(a3)") corr
            do i=1,n_reg1
            read (45 ,"(17x,f10.4)") q_mp2(i)
            end do
            else
            do i=1,n_reg1
            q_mp2(i) = q(i)
            end do
         end if
         open (44 , file="vdw.par")
      end if
      call readopt (iterld)
      iprint = 0
      call dg_ld (iterld,iprint)
      call solvout (iterld,do_gas)
      close (44 )
30    continue
      close (43 )
      close (45 )
1000  format(1x,A8,F8.1,2F10.4,3F9.4)
2000  format(i3)
      end
      Block Data
      implicit Real*8 (a-h, o-z)
      PARAMETER (MXATM=500 )
      parameter (mxcenter=50 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /born/ rg_reg1, rgim
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /pctimes/ ndxp,itl,itp
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /lra/ clgvn, slgvn
      data rg,drg,drg_inner,rdcutl,out_cut/26.0 ,3.0 ,1.0 ,6.1 ,16.0 /
      data dxp0/0.0 ,0.0 ,0.0 /
      data rgim/2.0 /
      data ndxp,itl,itp/10 ,399 ,5/
      data rzcut/1.1 /
      data rp/2.3 ,2.3 ,2.15 ,2.00 ,2.40 ,2.65 ,3.0 ,3.25 ,2.65 ,2.85 ,        3.2 ,3.0 ,2.32 ,2.65 ,2.8 ,2.46 ,2.5 ,2.58 ,1.82 ,1.70 ,        3.1 ,3.2 ,3.0 ,3.2 ,3.16 ,2.8 ,3.06 ,2.38 ,1.5 ,2.00 ,        1.91 ,1.89 ,1.93 ,1.84 ,1.57 ,1.50 ,1.88 ,1.55 ,2.00 ,2.50 ,        3.00 ,3.00 ,3.44 ,3.00 ,3.25 ,2.70 ,1.75 ,2.00 ,2.00 ,2.00 ,        2.00 ,2.00 ,2.00 ,2.00 ,2.25 ,1.98 ,2.45 ,2.44 ,3.00 ,3.70 ,        3.80 ,3.55 ,3.58 ,2.92 ,2.30 ,3.00 ,3.00 ,3.00 ,3.00 ,3.00 ,        3.00 ,1.77 ,1.94 ,2.00 ,2.55 ,3.00 ,3.00 ,3.00 ,3.00 ,3.00 ,        3.50 ,3.00 /
      end
      subroutine readopt (iterld)
      implicit Real*8 (a-h, o-z)
      PARAMETER (MXATM=500 )
      parameter (mxcenter=50 )
      character*1 dash(72 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /born/ rg_reg1, rgim
      common /aname/ atom(mxatm),molname,ssname
      common /pcgmcntr/ pcenter(3)
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /pctimes/ ndxp,itl,itp
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /lra/ clgvn, slgvn
      character*8 atom
      character*13 molname
      character*4 ssname
      character*2 rpinp
      data dash/72 *"-"/
      do i=1,82 
      vdwc6(i) = 0.0 
      end do
      read(44 ,1202 ) iterld, ndxp, dxp0(1), clgvn, slgvn,tds0
      read(44 ,1200 ) srp, rp(6), rp(7), rp(9), rp(10 )
      read(44 ,1200 ) rp(13 ), rp(14 ), rp(15 ), rp(22 ), rp(24 )
      read(44 ,1200 ) rp(16 ), rp(25 )
      read(44 ,1200 ) vdwc6(1), vdwc6(6), vdwc6(9), vdwc6(13 ),vdwc6(22 )
      read(44 ,1200 ) vdwc6(24 ), vdwc6(7)
      read(44 ,1200 ) vdwsl, phobsl, ephil1,ephil2, rzcut
      phobsl=phobsl/10.0d0 
      vdwsl=-vdwsl/10.0d0 
      do i=1, n_reg1
      if (vdwc6(iacw(i))  EQUAL 0.0d0 ) then
        if(vdwc6(iacw(i)-1) >< 0) then
        vdwc6(iacw(i)) = vdwc6(iacw(i)-1)
        else
        vdwc6(iacw(i)) = vdwc6(iacw(i)-2)
        end if
      end if
      if (vdwc6(iacw(i))  EQUAL 0.0d0 ) then
      if (iacw(i) <= 17 ) vdwc6(iacw(i))=vdwc6(6)
      if (iacw(i) > 17 ) vdwc6(iacw(i))=vdwc6(22 )
      end if
      end do
      do 60 i=1,n_reg1
      if(iacw(i) >< 1  AND iacw(i) >< 2) then
      rpi(i) = rp(iacw(i))
      else
      r2min = 100.0 
      jmin = 21 
      do 50 j=1,n_reg1
        if(j EQUAL i) go to 50 
        r2=(xw(1,i)-xw(1,j))**2+(xw(2,i)-xw(2,j))**2      +(xw(3,i)-xw(3,j))**2
        if(r2 < r2min) then
          r2min=r2
          jmin=j
        end if
50    continue
         if((iacw(jmin)) EQUAL 15 ) then
         iacw(i) = 2
         rpi(i) = rp(iacw(i))
         vdwc6(iacw(i)) = vdwc6(iacw(i)-1)
         else
         if(iacw(jmin) < 18 ) rpi(i) = srp * rp(iacw(jmin))
         if(iacw(jmin) >= 18 ) rpi(i) = (srp-0.1 )*rp(iacw(jmin))
         end if
      end if
60    continue
      read (45 ,"(a2)") rpinp
      if (rpinp EQUAL "rp") then
      read (45 ,*) nrp
      do j=1,nrp
      read (45 ,*) i, rpi(i)
      end do
      end if
      if (iprint EQUAL 1) then
      write(6,"(/,\"LD parameters: \",/)")
      write(6,1202 ) iterld, ndxp, dxp0(1), clgvn, slgvn
      write(6,1200 ) vdwsl, phobsl, ephil1, ephil2, rzcut
      end if
      write(6,95 )
      do i=1,n_reg1
      write (6,102 ) latom(i),(xw(j,i),j=1,3),q(i),q_gas(i),              iacw(i), rpi(i),vdwc6(iacw(i))
      enddo
      write(6,103 ) dash
 95   format(//,1x,"atom #",5x,"x",8x,"y",8x,"z",3x,"Q_pcm",2x,"Q_gas",1x,"atom_type",3x,"rp",4x,"VdWC6"/)
 102  format (i6,3f9.3,2f6.2,i6,7x,f4.2,4x,f4.2)
 103  format (/1x,72a1/)
1200  format (5f6.3)
1202  format (i2,i3,4f6.3)
      pcenter(1)=0.0 
      pcenter(2)=0.0 
      pcenter(3)=0.0 
      do i=1,n_reg1
         pcenter(1)=pcenter(1)+xw(1,i)/n_reg1
         pcenter(2)=pcenter(2)+xw(2,i)/n_reg1
         pcenter(3)=pcenter(3)+xw(3,i)/n_reg1
      enddo
      dmax=0.0d0 
      do i=1,n_reg1
      rg_reg1=(pcenter(1)-xw(1,i))**2 + (pcenter(2)-xw(2,i))**2 +        (pcenter(3)-xw(3,i))**2
      if (rg_reg1 > dmax) dmax=rg_reg1
      end do
      rg_reg1=sqrt(dmax)+2.4d0 
      rg=rg_reg1+rg
      rg_inner=rg_reg1+rgim+2.0d0 
      if (iprint EQUAL 1) then
      write(6,106 ) pcenter
      write(6,107 ) rg_reg1
      end if
      if (drg < 3.0 ) then
         write(6,996 )
         stop
      elseif (ndxp > mxcenter) then
         write(6,997 ) mxcenter
         stop
      endif
 106  format(/,1x,"Center of the cubic grid is defined as centroid of",     " solute atoms",//,20x,"center -",3f8.3)
 107  format(/,"Radius of the centroid: ",f6.2,//)
 996  FORMAT(//" PROGRAM EXIT: DRG SHOULD BE AT LEAST 3.O A"//)
 997  FORMAT(//" PROGRAM EXIT: MAXIMUN # OF GRIDS IS",I3//)
      return
      end
      subroutine dg_ld (iterld,iprint)
      implicit Real*8 (a-h,o-z)
      parameter (mxcenter=50 )
      PARAMETER (MXATM=500 )
      common /aname/ atom(mxatm),molname,ssname
      common /pctimes/ ndxp,itl,itp
      common /pcgmcntr/ pcenter(3)
      common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /atom_fs/ atomfs(mxatm)
      character*8 atom
      character*13 molname
      character*4 ssname
      integer i,j
      real*8 center_new(3)
      character*1 dash(72 )
      data dash/72 *"-"/
      esum = 0.0d0 
      evqdq = 0.0d0 
      ecor=0.0d0 
      do i=1,n_reg1
      atomfs(i)=0.0d0 
      end do
      do 10 i=1,ndxp
         call ran_shift(i,pcenter,center_new)
         call lgvnx(center_new,elgvn,eqm,ndipole,i,iterld,iprint)
         esum = esum + elgvn
         if (iterld EQUAL 1) then
         call sci_lgvn(elgwa,elgvni,tds,ndipole,i)
         temp_elgvn(i)=elgwa
         tdsl(i) = tds
         else
         elgwa= 0.0d0 
         temp_elgvn(i)=elgvn
         tdsl(i) = 0.0d0 
         end if
         tdsw_a = phobsl * fsurfa(i)
         if (iprint EQUAL 1) then
         write(6,"(/,\" Averaged atomic contributions to E_phobic:\",/)")
         do j=1,n_reg1
         write(6,112 ) j, atom(j), phobsl*atomfs(j)/float(i)
         end do
         end if
         if (iterld EQUAL 1) write(6,105 ) elgwa,evdwl(i),tdsw_a,-tds+tdsw_a
         if (iterld EQUAL 0) write(6,105 ) elgvn,evdwl(i),tdsw_a,-tds+tdsw_a
         write(6,102 ) dash
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
 102  format(1x,72a1)
 105  format(/1x,"Elgvn = ",f8.3,"   Evdw = ",f8.3,       "   -TdS_phobic = ",f8.3,"   -TdS_total = ",f8.3,/)
 112  format(i3,5x,a8,f9.3)
      end
      subroutine ran_shift(i,center1,center2)
      implicit Real*8 (a-h,o-z)
      parameter (mxcenter=50 )
      PARAMETER (MXATM=500 )
      common /pctimes/ ndxp,itl,itp
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
      dimension oshift(3*mxcenter)
      dimension dxp(3), center2(3), center1(3)
      save oshift
       if (i EQUAL 1) then
         iseed = -931 
         idum = 1
         dumm = ran2(iseed)
         do kk = 1, 3*ndxp
         oshift(kk) = ran2 (idum)
         end do
       end if
       fact=drg_inner
       if(rg_inner EQUAL 0.0d0 ) fact=drg
       if(i EQUAL 1)fact=0.0d0 
       dxp(1)=fact*(1.0d0 -2.0d0 *oshift(3*i-2))
       dxp(2)=fact*(1.0d0 -2.0d0 *oshift(3*i-1))
       dxp(3)=fact*(1.0d0 -2.0d0 *oshift(3*i))
       if(i >< 1)then
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
       write(6,100 ) center2
       return
100    format(/ " Grid origin            ",3f9.3)
       end
      subroutine lgvnx(center1,elgvn,eqm,ndipole,ientro,iterld,iprint)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000 )
      parameter (mxatm=500 )
      parameter (mxcenter=50 )
      PARAMETER (ONE=1.0d0 , ZERO=0.0d0 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common/scrat5/  rz1(mxlgvn),rz_vdw(mxlgvn), iz(mxlgvn)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /atom_fs/ atomfs(mxatm)
      common /outer_surf/ isd (mxlgvn)
      common /lra/ clgvn, slgvn
      real*8 center1(3)
      integer*2 isd
      real*8 elgvn,elgvna
      real*8 fs, gri_sp, efn, fma
      real*8 vdwsur(mxatm)
      character*1 dash(72 )
      data dash/72 *"-"/
      elgvn=0.0d0 
      sres= 0.0d0 
      elgvn=0.0d0 
      fsurfa(ientro)=0.0d0 
      evdwl(ientro)=0.0d0 
      ndipole=0
      idum = 1
      tds = 0.0d0 
      call gen_gridx(center1,ndipole,ientro,0,iprint)
      do i=1,ndipole
      xd(1,i) = xl(1,i)
      xd(2,i) = xl(2,i)
      xd(3,i) = xl(3,i)
      end do
      call ef_ld(ndipole,1)
      efn_max=-10.0d0 
      elgvn = 0.0d0 
      do i = 1,n_reg1
      vdwsur(i) = 0.0d0 
      end do
      do 10 i=1,ndipole
      gri_sp=drg_inner
      if(i > n_inner) gri_sp=drg
      if(i EQUAL n_inner+1) elgvn_i = elgvn
      efn=dsqrt(da(1,i)*da(1,i)+da(2,i)*da(2,i)+da(3,i)*da(3,i))
      if (efn > efn_max) efn_max=efn
      call vlgvn(efn,elgvn,xjunk,fma,gri_sp)
      xmua(1,i)=fma*da(1,i)/efn
      xmua(2,i)=fma*da(2,i)/efn
      xmua(3,i)=fma*da(3,i)/efn
      if (i > n_inner) goto 10 
      if (rz_vdw(i) <= rzcut) then
        epot = 0.0d0 
        do k=1,n_reg1
        rx = xl(1,i) - xw(1,k)
        ry = xl(2,i) - xw(2,k)
        rz = xl(3,i) - xw(3,k)
        rqd = sqrt (rx*rx+ry*ry+rz*rz)
        epot = epot + q(k)/rqd
        end do
       if (epot < 0.0d0 ) epot=-epot
        if (epot <= ephil1) then
         fs=1.0d0 
        end if
        if((epot > ephil1) AND (epot <= ephil2)) then
         fs=1.0d0 -((epot-ephil1)/(ephil2-ephil1))*(1.0 -sres)
         elseif(epot > ephil2) then
         fs=sres
        endif
        atomfs(iz(i)) = atomfs(iz(i))+fs
        fsurfa(ientro)=fsurfa(ientro)+fs
        vdwsur(iz(i)) = vdwsur(iz(i)) + 1.0d0 
      endif
      ! #GOTO 10 
 10   continue
      evdwl(ientro) = 0.0d0 
      do k=1,n_reg1
      evdwl(ientro) = evdwl(ientro) + vdwc6(iacw(k))*vdwsur(k)
      end do
      evdwl(ientro) = vdwsl*evdwl(ientro)
      elgvn = clgvn * elgvn
      elgvn_i = clgvn * elgvn_i
      write(6,1001 ) ndipole,elgvn, elgvn_i
      if (iterld. eq. 0) return
      call ef_ld(ndipole,0)
      elgvna=0.0 
      do 20 i=1,ndipole
      gri_sp=drg_inner
      if(i > n_inner) gri_sp=drg
      efn=dsqrt(da(1,i)*da(1,i)+da(2,i)*da(2,i)+da(3,i)*da(3,i))
      call vlgvn(efn,elgvna,dumm,fma,gri_sp)
      if(i <= n_inner  OR  isd (i) EQUAL 1  OR  mod(i,3) EQUAL 0) then
      ddd=3.0d0 /(sqrt(rz1(i))+2.0d0 )
      xmua(1,i)=ddd*fma*da(1,i)/efn
      xmua(2,i)=ddd*fma*da(2,i)/efn
      xmua(3,i)=ddd*fma*da(3,i)/efn
      go to 20 
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
 1001 format(" Total number of Langevin dipoles      : ",i10/       " Noniterative lgvn energy (total,inner): ",2f10.3/)
      end
      function ran2 (idum)
      implicit Real*8 (a-h,o-z)
      Parameter (m=714025 , ia=1366 , ic=150889 , rm=1.0 /m)
      Integer ir(97 )
      save ir, iy
      Data iff /0/
      if (idum < 0 OR iff EQUAL 0) then
        iff = 1
        idum=mod(ic-idum,m)
        do 11 j=1,97 
          idum=mod(ia*idum+ic,m)
          ir(j) = idum
11      continue
        idum=mod(ia*idum+ic,m)
        iy=idum
      end if
      j=1+(97 *iy)/m
      if(j > 97 OR j < 1) then
      print*,"Problems with the random number generator, j = ", j
      stop
      end if
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j) = idum
      return
      end
      subroutine ef_ld (ndipole,idiel)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000 )
      parameter (mxatm=500 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      if (idiel EQUAL 0) then
      do 5 i=1,3
      do 5 j=1,ndipole
5     da(i,j) = 0.0d0 
      do 10 j=1,ndipole
      do 10 i=1,n_reg1
      ri = xl(1,j)-xw(1,i)
      rj = xl(2,j)-xw(2,i)
      rk = xl(3,j)-xw(3,i)
      r2 = ri*ri + rj*rj + rk*rk
      r3 = r2 * dsqrt(r2)
      qr = q(i) / r3
      da(1,j) = da(1,j) + qr*ri
      da(2,j) = da(2,j) + qr*rj
      da(3,j) = da(3,j) + qr*rk
10    continue
      end if
      if (idiel EQUAL 1) then
      do 15 i=1,3
      do 15 j=1,ndipole
15    da(i,j) = 0.0d0 
      do 20 j=1,ndipole
      do 20 i=1,n_reg1
      ri = xl(1,j)-xw(1,i)
      rj = xl(2,j)-xw(2,i)
      rk = xl(3,j)-xw(3,i)
      r2 = ri*ri + rj*rj + rk*rk
      r1 = dsqrt(r2)
      r3 = r2 * r1
      ddd = 1.7d0 /sqrt(r1+2.0 )
      qr = ddd * q(i) / r3
      da(1,j) = da(1,j) + qr*ri
      da(2,j) = da(2,j) + qr*rj
      da(3,j) = da(3,j) + qr*rk
20    continue
      end if
      return
      end
      subroutine sci_lgvn(energy,elgvni,tds,nd,icent)
      implicit Real*8 (a-h,o-z)
      parameter (mxatm=500 )
      parameter (mxlgvn=10000 )
      parameter (mxcenter=50 )
      parameter (mxpair2=10000000 )
      parameter (mxpair=5000000 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /pcefaefb/efa(3,mxlgvn)
      common efal (3,mxlgvn)
      common /pcdipcut/ rdcutl,out_cut
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      common /pctimes/ ndxp,itl,itp
      integer*2 jp, jp2
      common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn),                 jp2(mxpair2),ip3(0:mxlgvn)
      common /lra/ clgvn, slgvn
      real*8 energy
      integer nd
      integer i,l,icent,iopen
      real*8 conv2
      real*8 epom(44 )
      character*1 dash(72 )
      data dash/72 *"-"/
      data iopen/1/
      tds = 0.0d0 
      write(6,1006 )
      call pairlistw(nd,xd)
      write(6,1008 )
      conv2=-332.0d0 *slgvn
      do i=1,10 
      epom(i) = 0.0d0 
      end do
      do 40 l=1,itl
      if ( %FALSE% ) then
      open (50 , file="cs_dipoles.xyz")
      nd2plot=0
      do i=1,nd
         if (xl(1,i) EQUAL 0.0 ) nd2plot=nd2plot+1
      enddo
      write(50 ,2006 ) nd2plot
      write(50 ,*)
      do i=1,nd
         scale = 25.0 
         if (i <= n_inner) scale = 50.0 
         if (xl(1,i) EQUAL 0.0 ) then
         write (50 ,2007 ) xl(1,i),xl(2,i), xl(3,i),   scale*xmua(1,i),scale*xmua(2,i),scale*xmua(3,i)
         endif
      enddo
      end if
 2006 format(i5)
 2007 format("X",6f10.3)
      eita=0.0 
      do i=1,n_inner
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eiti=conv2*eita
      do i=n_inner+1,nd
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eita=conv2*eita
      write(6,1013 ) l-1,eita,eiti,eita-tds
      if (l EQUAL 1) elgvni=eita
      ni=10 
      if(l > 5) then
        eave=0.0d0 
        epom(ni+1)= eita - tds
        do i=1,ni
        epom(i) = epom(i+1)
        eave=eave+epom(i)
        end do
        eave=eave/float(ni)
        iconvp=1
        iconv=1
        do i=1,ni
        if((100.0d0 *abs((eave-epom(i))/eave)) > 0.1d0 ) iconv=0
        iconvp = iconvp * iconv
        end do
        if(l > 16 AND epom(ni) > epom(ni-1) AND   epom(ni-1) > epom(ni-2))  iconvp = 1
        if(iconvp EQUAL 1) then
        goto 50 
        end if
      end if
      call newf_lcut(nd,l)
      if(l < 10 OR l EQUAL 45 ) call updatelong(nd,l)
      stepaw=0.2 
      if(l < 4) stepaw=0.10 
      if(l < 11 ) stepaw=0.20 
      if(l > 28 ) stepaw=0.30 
      if(l > 80 ) stepaw=0.50 
      call mu_mu_l (nd,stepaw,tds)
 40   continue
      eita=0.0d0 
      do i=1,n_inner
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eiti=conv2*eita
      do i=n_inner+1,nd
      eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
      enddo
      eita=conv2*eita
      write(6,1014 ) l-1,eita,eiti,eita-tds
 50   continue
      energy=eita
      return
 1002 format(1x,"Total number of dipoles = ",i8)
 1006 format(/,1x,"Iterative LD calculation for the chosen grid")
 1008 format(//1x,"Iteration",9x,"Elgvn",5x,"Elgvn_inner",5x,      "Elgvn-TdS_immob"/)
 1013 format(5x,i3,4x,f12.3,5x,f12.3,5x,f12.3)
 1014 format(5x,i3,4x,f12.3,5x,f12.3,5x,f12.3//)
      end
      subroutine vlgvn (efn,ea,tds,fma,gri_sp)
      implicit Real*8 (a-h,o-z)
      common /lra/ clgvn, slgvn
      xdrg= (gri_sp/3.0d0 )**3
      conv=-332.0d0 *slgvn
      dipmax=0.29d0 
      dip=sqrt(xdrg)*dipmax
      aktm=332.0d0 /0.6d0 
      dddi=3.0 
      x=dipmax*efn*aktm/dddi
      x2=x*x
      x4=x2*x2
      x6=x2*x4
      ex=exp(x)
      exm=1.0d0 /ex
      algvn=(ex+exm)/(ex-exm) - 1.0d0 /x
      fma=dip*algvn
      tds = -dlog(3.14159d0 /(2.0d0 *acos(algvn)))
      corrf = 5d0 *atan(x/27.0d0 )*exp(-x6/32.0d0 )
      tds = tds + corrf
      tds = tds * 9.0d0 *xdrg**(2.0d0 /3.0d0 )
      tds = tds/17.3d0 
      ea=ea+conv*fma*efn
      return
      end
      subroutine newf_lcut(nd,l)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000 )
      parameter (mxpair2=10000000 )
      parameter (mxpair=5000000 )
      common /pcefaefb/efa(3,mxlgvn)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      integer*2 jp, jp2
      common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn),                 jp2(mxpair2),ip3(0:mxlgvn)
      integer nd,l
      integer i,i1,index,i2
      real*8 d2
      real*8 c5,rmu3a,rmu3ar
      integer npair,np
      real*8 r1,r2,r3
      do i=1,nd
         efa(1,i)=da(1,i)
         efa(2,i)=da(2,i)
         efa(3,i)=da(3,i)
      enddo
      npair=0
      do i1=1,nd-1
         np=ip(i1)-ip(i1-1)
         if(np > 0) then
            do index=1,np
               npair=npair+1
               i2=jp(npair)
               r1=xd(1,i1)-xd(1,i2)
               r2=xd(2,i1)-xd(2,i2)
               r3=xd(3,i1)-xd(3,i2)
               d2=r1*r1+r2*r2+r3*r3
               c5=-1.0 /(d2*d2*dsqrt(d2))
               rmu3a=3.0 *(r1*xmua(1,i2)+r2*xmua(2,i2)+r3*xmua(3,i2))
               rmu3ar=3.0 *(r1*xmua(1,i1)+r2*xmua(2,i1)+r3*xmua(3,i1))
               efa(1,i1)=efa(1,i1)+c5*(xmua(1,i2)*d2-rmu3a*r1)
               efa(2,i1)=efa(2,i1)+c5*(xmua(2,i2)*d2-rmu3a*r2)
               efa(3,i1)=efa(3,i1)+c5*(xmua(3,i2)*d2-rmu3a*r3)
               efa(1,i2)=efa(1,i2)+c5*(xmua(1,i1)*d2-rmu3ar*r1)
               efa(2,i2)=efa(2,i2)+c5*(xmua(2,i1)*d2-rmu3ar*r2)
               efa(3,i2)=efa(3,i2)+c5*(xmua(3,i1)*d2-rmu3ar*r3)
            enddo
         endif
      enddo
         return
         end
      subroutine updatelong(nd,l)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000 )
      parameter (mxpair2=10000000 )
      parameter (mxpair=5000000 )
      common efal(3,mxlgvn)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      integer*2 jp, jp2
      common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn),                 jp2(mxpair2),ip3(0:mxlgvn)
      integer nd
      integer i,l,npair,np,index,j
      real*8 d2,r1,r2,r3
      real*8 c5,rmu3a,rmu3ar
      do j=1,3
         do i=1,nd
            efal(j,i)=0.0d0 
         enddo
      enddo
      npair=0
      do i=1,nd-1
         np=ip2(i)-ip2(i-1)
         if(np > 0) then
            do index=1,np
               npair=npair+1
               j=jp2(npair)
               r1=xd(1,i)-xd(1,j)
               r2=xd(2,i)-xd(2,j)
               r3=xd(3,i)-xd(3,j)
               d2=r1*r1+r2*r2+r3*r3
               c5=-1.0d0 /(d2*d2*dsqrt(d2))
               rmu3a=3.0 *(r1*xmua(1,j)+r2*xmua(2,j)+r3*xmua(3,j))
               rmu3ar=3.0 *(r1*xmua(1,i)+r2*xmua(2,i)+r3*xmua(3,i))
               efal(1,i)=efal(1,i)+c5*(xmua(1,j)*d2-rmu3a*r1)
               efal(2,i)=efal(2,i)+c5*(xmua(2,j)*d2-rmu3a*r2)
               efal(3,i)=efal(3,i)+c5*(xmua(3,j)*d2-rmu3a*r3)
               efal(1,j)=efal(1,j)+c5*(xmua(1,i)*d2-rmu3ar*r1)
               efal(2,j)=efal(2,j)+c5*(xmua(2,i)*d2-rmu3ar*r2)
               efal(3,j)=efal(3,j)+c5*(xmua(3,i)*d2-rmu3ar*r3)
            enddo
         endif
      enddo
      return
      end
      subroutine mu_mu_l (n,step_,tds)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000 )
      parameter (mxatm = 500 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /pcefaefb/efa(3,mxlgvn)
      common efal(3,mxlgvn)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /outer_surf/ isd (mxlgvn)
      integer*2 isd
      tds = 0.0d0 
      tds_sum=0.0d0 
      do 10 i=1,n
          if(isd(i) EQUAL 0) then
          efa(1,i)=efa(1,i)+efal(1,i)
          efa(2,i)=efa(2,i)+efal(2,i)
          efa(3,i)=efa(3,i)+efal(3,i)
          efna=efa(1,i)*efa(1,i)+efa(2,i)*efa(2,i)+efa(3,i)*efa(3,i)
          efna=dsqrt(efna)
          gri_sp=drg_inner
          if(i > n_inner) gri_sp=drg
          call vlgvn (efna,xjunka,tds,fma,gri_sp)
          tds_sum = tds_sum + tds
          xmua(1,i)=xmua(1,i)+step_*(fma*efa(1,i)/efna-xmua(1,i))
          xmua(2,i)=xmua(2,i)+step_*(fma*efa(2,i)/efna-xmua(2,i))
          xmua(3,i)=xmua(3,i)+step_*(fma*efa(3,i)/efna-xmua(3,i))
          end if
10    continue
      tds=tds_sum
      return
      end
      subroutine pairlistw(nd,dipx)
      implicit Real*8 (a-h,o-z)
      parameter (mxatm=500 )
      parameter (mxpair=5000000 )
      parameter (mxpair2=10000000 )
      parameter (mxlgvn=10000 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn),                 jp2(mxpair2),ip3(0:mxlgvn)
      common /outer_surf/ isd (mxlgvn)
      integer*2 jp,jp2,isd
      integer   ip,ip2,ip3
      integer nd
      real*8 dipx(3,mxlgvn)
      integer i,j
      integer npair,npair2,npair3,ntotl
      real*8 d2,rdcut2,r1,r2,r3,out_cut2
      integer mask1
      real*8 rcoff,rcoff2,out_cut
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
            if(isd(i) EQUAL 1.0 and. isd(j) EQUAL 1) goto 15 
            r1=dipx(1,i)-dipx(1,j)
            r2=dipx(2,i)-dipx(2,j)
            r3=dipx(3,i)-dipx(3,j)
            d2=r1*r1+r2*r2+r3*r3
            rcoff=2.5 
            rcoff2=rcoff*rcoff
            mask1=d2/rcoff2
            if(mask1 >= 1) then
               ntotl=ntotl+1
               if(d2 < rdcut2) then
                  npair=npair+1
                  if (npair > mxpair) then
                     write(6,1000 )npair,mxpair
                     stop
                  endif
                  jp(npair)=j
               else if(d2 >= rdcut2 AND d2 < out_cut2)then
                  npair2=npair2+1
                  if (npair2 > mxpair2) then
                     write(6,1002 ) npair2,mxpair2
                     stop
                  endif
                  jp2(npair2)=j
               endif
            endif
 15      continue
         ip(i)=npair
         ip2(i)=npair2
         ip3(i)=npair3
 10   continue
      write(6,1001 ) ntotl,npair,rdcutl,npair2,out_cut,            npair3,out_cut
      return
 1000 format(/,1x,i8," exceeds program mxpair dimension - ",i8)
 1002 format(/,1x,i8," exceeds program mxpair2 dimension - ",i8)
 1001 format(/1x,"dipole pair list"/     /,1x,"       total possible pairs  - ",i10," out of ",          "   2.500 A cutoff",     /,1x,"       pairs selected        - ",i10," within ",          f8.3," A cutoff",     /,1x,"       pairs selected (long) - ",i10," within ",          f8.3," A cutoff",     /,1x,"       pairs selected (long) - ",i10," out of ",          f8.3," A cutoff")
      end
      subroutine solvout (iterld, do_gas)
      implicit Real*8 (a-h,o-z)
      PARAMETER (MXATM=500 )
      parameter (mxcenter=50 )
      common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /aname/ atom(mxatm),molname,ssname
      character*8 atom
      character*13 molname
      character*4 ssname
      logical do_gas
      character*1 dash(72 )
      data dash/72 *"-"/
      write(6,20 ) dash
      erelax = evqdq
      dgsolv=elgwa+etds +evdw+ebw
      if (iterld EQUAL 0) dgsolv = 0.0d0 
      dgsolvni=elgvn+ephob+evdw+ebw
      write(6,100 )
      if (iterld EQUAL 0) write(6,85 )
      if (iterld EQUAL 1) write(6,101 ) elgwa
      if (iterld EQUAL 0) write(6,101 ) elgvn
      write(6,102 ) etds
      write(6,103 ) evdw
      write(6,104 ) ebw
      if (do_gas) write(6,105 ) erelax
      if (iterld EQUAL 1) write(6,106 ) dgsolv
      if (iterld EQUAL 0) write(6,106 ) dgsolvni
      if (iprint EQUAL 1) then
      write(6,110 )
      write(6,115 ) elgvn
      write(6,117 ) dgsolvni
      end if
      if (iterld EQUAL 1) then
      open (43 ,file="cs.arc",access="append")
      if (do_gas) then
      write(6,199 )
      write(6,200 ) molname,ssname,elgwa, evdw, etds , erelax, ebw,             dgsolv-etds,dgsolv
      write(43 ,201 ) molname,ssname,elgwa, evdw, etds , erelax, ebw,             dgsolv-etds,dgsolv
      else
      write(6,299 )
      write(6,300 ) molname,ssname,elgwa, evdw, etds , ebw,             dgsolv-etds,dgsolv
      write(43 ,301 ) molname,ssname,elgwa, evdw, etds , ebw,             dgsolv-etds,dgsolv
      end if
      close(43 )
      end if
199   format(//,"Molecule               lgvn    VdW   -TdS   Relax",          "    Born   dHsolv   dGsolv")
200   format (a13,1x,a4,f9.1,3f7.1,3f9.1/)
201   format (a13,1x,a4,f9.1,3f7.1,3f9.1)
299   format(//,"Molecule               lgvn    VdW   -TdS        ",          "    Born   dHsolv   dGsolv")
300   format (a13,1x,a4,f9.1,2f7.1,7x,3f9.1/)
301   format (a13,1x,a4,f9.1,2f7.1,7x,3f9.1)
      return
 20   FORMAT(1X,72A1//" FINAL RESULTS")
  85  format(" Noniterative LD results:"/,"Warning! Noniterative LD is not parametrized in cs2.1"/,"         Put iterld=1 in the vdw.par file to run iterative LD"/)
 100   format(//," *****************************"/           " Langevin Dipoles Calculation "/           " *****************************"//  " Transfer of solute at 298K: GAS,1M -> WATER,1M (kcal/mol)",/)
 101   format(" Langevin energy         ",10x,f20.2)
 102   format(" -TdS                    ",10x,f20.2)
 103   format(" VdW energy              ",10x,f20.2)
 104   format(" dBorn+dOnsager          ",10x,f20.2)
 105   format(" Solute relaxation       ",10x,f20.2)
 107   format(" Correlation correction  ",10x,f20.2)
 106   format(" Total dG                ",10x,f20.2/)
 110   format(//," ****************************"/           " Noniterative LD calculation "/           " ****************************"//       " Transfer of a solute at 298K: GAS,1M -> WATER,1M",//)
 115   format(" Noniterative lgvn energy",10x,f20.2)
 117   format(" dG                      ",10x,f20.2)
       end
      subroutine vbornx(nd_lgvn,eb,center1)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=5000 )
      parameter (mxatm=500 )
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /born/ rg_reg1, rgim
      character*1 dash(72 )
      dimension center1(3)
      data dash/72 *"-"/
       write(6,201 ) dash
       chqm = 0.0d0 
       dqmx = 0.0d0 
       dqmy = 0.0d0 
       dqmz = 0.0d0 
       conv     =  -(1.0d0 -1.0d0 /80.0d0 )
       conv_dip =  -158.0d0 /161.0d0 
       v_outer = 27.0 * float(nd_lgvn - n_inner)
       b=rg_reg1+rgim
       b=(3.0 *v_outer*0.25 /3.14159 + b*b*b)**(1.0 /3.0 )
       b3=b*b*b
       do 10 i=1,n_reg1
       chqm = chqm + q(i)
       dqmx=dqmx+q(i)*(xw(1,i)-center1(1))
       dqmy=dqmy+q(i)*(xw(2,i)-center1(2))
       dqmz=dqmz+q(i)*(xw(3,i)-center1(3))
10     continue
       dqm = dsqrt(dqmx**2 + dqmy**2 + dqmz**2)
       if(abs(chqm) > 1.02 ) fch = 1.0d0 
       if(abs(chqm) > 0.02 ) fch = 0.1d0 
       if(abs(chqm) <= 0.02 ) fch = 0.90d0 
       if(dqm >= 5.5d0 ) fch=0.0d0 
       eda = 166.0d0 *conv_dip*((fch*dqm)**2)/b3
       dqm = 4.8023 *dqm
       eb=eda
       eba = 0.0d0 
       if(abs(chqm) < 1.98 ) eba=166.0d0 *conv*chqm**2/b
       if(abs(chqm) >= 1.98 ) eba=166.0d0 *conv*0.75d0 *chqm**2/b
       eb=eb+eba
       write (6,200 ) b,dqm,chqm,eba,eda,eb
       return
 200   format(//" Born radius equals - ",f9.2,1x,"A"//            " molecular dipole (Debye)        - ",f9.2/            " molecular charge                - ",f9.2/            " Born  monopole energy           - ",f9.2/            " scaled Onsager dipole energy    - ",f9.2//            " total bulk energy difference    - ",f9.2/)
 201   format(/1x,72a1//," Estimation of the bulk energy using Born ",            "equation")
       end
      subroutine vatom (cor,vqdq,ndipole,iterld,iprint)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000 )
      parameter (mxatm=500 )
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      common /lra/ clgvn, slgvn
      character*1 dash(72 )
      dimension vq(mxatm)
      dimension rl(3)
      data dash/72 *"-"/
      if (iprint EQUAL 1) then
      write(6,201 ) dash
      write(6,202 )
      end if
      vqq = 0.0d0 
      vqdq = 0.0d0 
      cor = 0.0 
      do 20 j=1,n_reg1
          vq(j)=0.0d0 
          do 10 i=1,ndipole
             r2=0.0d0 
             sp=0.0d0 
             do 101 k=1,3
                rl(k)=-xl(k,i)+xw(k,j)
                sp=sp+rl(k)*xmua(k,i)
                r2=r2+rl(k)**2
 101         continue
             r1=dsqrt(r2)
             if (iterld EQUAL 1) then
             vq(j)=vq(j)+332 *sp/(r2*r1)
             else
             vq(j)=vq(j)+332 *(1.7d0 /dsqrt(r1+2.0D0 ))*sp/(r2*r1)
             end if
 10       continue
          dq_gas = q_gas(j) - q(j)
          dq_mp2 = q_mp2(j) - q(j)
          if (iprint EQUAL 1) write(6,300 ) j, zan(j), q(j), vq(j),                     vq(j)*q(j), vq(j)*dq_gas, vq(j)*dq_mp2
          vqq = vqq + vq(j)*q(j)
          vqdq = vqdq + vq(j)*dq_gas
          cor = cor + vq(j)*dq_mp2
 20   continue
      vqq = slgvn*vqq
      if(iterld EQUAL 0) vqq = clgvn*vqq
      vqdq = 0.5d0 *vqdq
      if(iterld EQUAL 0) vqdq = 0.8 *vqdq
      if (iprint EQUAL 1) then
      write(6,400 ) vqq
      write(6,500 ) vqdq
      write(6,600 ) cor
      end if
 201  format(1x,72a1/," Potentials at solute atoms [kcal/mol]",/)
 202  format(2x,"Atom #",2x,"Type",2x,"Charge",6x,"Vq     ",       2x,"Vq*Q",2x,"Vq*dQ_pcm",2x,"Vq*dQ_mp2")
 300  format(2x,i5,3x,f4.1,2x,f5.2,1x,4(2x,f7.1))
 400  format(/,"Langevin energy calculated from potential at solute",      " atoms:",3x,f9.2,2x,"kcal/mol")
 500  format("Solute relaxation energy calculated from PCM charges:",      "       ",f9.2,2x,"kcal/mol")
 600  format("Correction for electron correlation effects (dG_cor):",      "       ",f9.2,2x,"kcal/mol")
      end
      subroutine elgvn_ave(iterld,ncenter)
      implicit Real*8 (a-h,o-z)
      parameter (mxcenter=50 )
      common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /volume/ nvol(mxcenter)
      character*1 dash(72 )
      data dash/72 *"-"/
      eav=0.0d0 
      eboltz=0.0d0 
      enorm=0.0d0 
      write(6,1000 )
      temper=300.0d0 
      wlowest=10000.0d0 
      do 10 jj=1,ncenter
      if(temp_elgvn(jj) < wlowest) then
         wlowest=temp_elgvn(jj)
         nlowest=jj
      endif
 10   continue
      elowest=temp_elgvn(nlowest)
      do 20 i=1,ncenter
         write(6,1001 ) i,temp_elgvn(i)
         eav=eav+temp_elgvn(i)/dble(ncenter)
         expe=dexp(-(temp_elgvn(i)-elowest)/0.002d0 /temper)
         eboltz=eboltz+temp_elgvn(i)*expe
         enorm=enorm+expe
 20      continue
      eboltz=eboltz/enorm
      write(6,1002 )elowest,eav,eboltz,temper
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
      write(6,"(/,\" Contributions to -TdS (kcal/mol)\")" )
      write(6,"(' Change in free-volume:  ',5x,f10.3)")  tds0
      write(6,"(' Hydrophobic          :  ',5x,f10.3)")  ephob
      write(6,"(' Dipolar saturation   :  ',5x,f10.3)")  -etds
      etds = tds0 + ephob - etds
      if (iprint EQUAL 1) then
      write(6,"(/,'Average VDW energy:            ',f10.3)") evdw
      write(6,"('Average -TdS (kcal/mol):       ',f10.3)") etds
      write(6,"('Average solvation entropy (eu):',f10.3)")          -1000.0d0 *etds/temper
      write(6,"('Average relaxation energy:     ',f10.3)") evqdq
      end if
      if (iterld EQUAL 1) elgwa = eav
 1000  format(1x," Final LD energies for different grids"//        "                    Elgvn     "/)
 1001  format(1x,i10,2x,f10.3)
 1002  format(//" lowest energy        --->  ",f10.3/          " mean energy          --->  ",f10.3/          " boltzmann <energies> --->  ",f10.3/          " for temperature      --->  ",f10.1," K"/)
       return
       end
      Real*8 function entropy (m)
      implicit real*8 (a-h,o-z)
      real*8 NA, kB, m, n
      NA = 6.023d23 
      kB = 1.381d%23 
      h = 6.626d%34 
      pi = 3.14159 
      T = 298.15 
      V = 1d%3 
      e = 2.71828 
      amu = 1.66d%27 
      n = 1.0 
      m = amu*m
      NA = NA*n
      S1 = sqrt((2*pi*m*kB*T)**3)/h**3
      S2 =  V*sqrt(e**5)/NA
      S = S1*S2
      S = dlog(S)
      S = NA*kB*S/n
      entropy = T*S/(4.18d0 *1000.0d0 )
      return
      end
      subroutine gen_gridx (center1,nld,ientro,iflag,iprint)
      implicit Real*8 (a-h,o-z)
      parameter (mxatm=500 )
      parameter (mxlgvn=10000 )
      parameter (mxcenter=50 )
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82 ),vdwc6(82 ),            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)             ,q_gas(mxatm),q_mp2(mxatm)
      common /born/ rg_reg1, rgim
      common /pcgrid/ drg,rg,dxp0(3),              rg_inner,drg_inner
      common /pcgmcntr/ pcenter(3)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),                da(3,mxlgvn)
      common/scrat5/  rz1(mxlgvn),rz_vdw(mxlgvn), iz(mxlgvn)
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /aname/ atom(mxatm),molname,ssname
      common /outer_surf/ isd (mxlgvn)
      common /volume/ nvol(mxcenter)
      character*8 atom
      character*13 molname
      character*4 ssname
      integer i0,i1,nld,iflag
      integer*2 isd
      real*8 center1(3)
      integer i,index,ii,jj,kk,no,n_out1
      integer ipro,jpro,kpro,igrid,jgrid,kgrid
      integer icube(0:133,0:133,0:133)
      real*8 rg2,rgi2,rgim2,drg2,drgi2,drg_i,drg_inner_i
      real*8 rp2(mxatm)
      real*8 ri,rj,rk,xloca,yloca,zloca,vdw_in,vdw_out
      real*8 xp(3),d2,xdrg,r6,d6,r3,d3
      nld=0
      ns=0
      i0=1
      i1=n_reg1
      rg2=rg*rg
      drg2=drg*drg
      drgi2=drg_inner*drg_inner
      rgi2=rg_inner*rg_inner
      rgim2=rgim*rgim
      drg_inner_i=1.0d0 /drg_inner
      drg_i=1.0d0 /drg
      do i=1,n_reg1
      rp2(i)=rpi(i)*rpi(i)
      end do
      limit_inner=int(rg_inner*2/drg_inner)
      limit_outer=int(rg*2/drg)
      if(mod(limit_inner,2) EQUAL 0) limit_inner=limit_inner+1
      if(mod(limit_outer,2) EQUAL 0) limit_outer=limit_outer+1
      mid_inner=int(limit_inner/2.0d0 +0.5d0 )
      mid_outer=int(limit_outer/2.0d0 +0.5d0 )
      if(rg_inner > 0.0d0 ) then ! build inner grid
      do 10 ii=0,limit_inner
         do 10 jj=0,limit_inner
            do 10 kk=0,limit_inner
               icube(ii,jj,kk)=0
               ri=ii-mid_inner
               rj=jj-mid_inner
               rk=kk-mid_inner
               d2=(ri*ri+rj*rj+rk*rk)*drgi2 !radius from center
               if(d2 <= rgi2)icube(ii,jj,kk)=1   !trim to a sphere
 10   continue
      nvol(ientro) = 0
      do 20 i=i0,i1 !loop over solute atoms
         ipro=int(mid_inner+(xw(1,i)-center1(1))*drg_inner_i)
         jpro=int(mid_inner+(xw(2,i)-center1(2))*drg_inner_i)
         kpro=int(mid_inner+(xw(3,i)-center1(3))*drg_inner_i)
         if(ipro >= 0 AND ipro <= limit_inner AND       jpro >= 0 AND jpro <= limit_inner AND       kpro >= 0 AND kpro <= limit_inner) then !within grid/dimensi
            do 30 ii=1,9
               igrid=ipro+ii-5
               do 30 jj=1,9
                  jgrid=jpro+jj-5
                  do 30 kk=1,9
                     kgrid=kpro+kk-5
                     if(icube(igrid,jgrid,kgrid) EQUAL 1)then
                        xloca=(igrid-mid_inner)*drg_inner+center1(1)
                        yloca=(jgrid-mid_inner)*drg_inner+center1(2)
                        zloca=(kgrid-mid_inner)*drg_inner+center1(3)
                        ri=xloca-xw(1,i)
                        rj=yloca-xw(2,i)
                        rk=zloca-xw(3,i)
                        d2=ri*ri+rj*rj+rk*rk !distance from grid point
                        if(d2 < rp2(i)) then
                           icube(igrid,jgrid,kgrid)=0
                           nvol(ientro) = nvol(ientro) + 1
                        end if
                     endif
 30         continue
         endif
 20   continue
      do 40 ii=0,limit_inner
         do 40 jj=0,limit_inner
            do 40 kk=0,limit_inner
               if(icube(ii,jj,kk) >< 0) then
                  xp(1)=center1(1)+(ii-mid_inner)*drg_inner
                  xp(2)=center1(2)+(jj-mid_inner)*drg_inner
                  xp(3)=center1(3)+(kk-mid_inner)*drg_inner
                  d2_min=10000.0d0 
                  do 60 i=i0,i1
                     ri=xp(1)-xw(1,i)
                     rj=xp(2)-xw(2,i)
                     rk=xp(3)-xw(3,i)
                     d2=ri*ri+rj*rj+rk*rk
                     d2=sqrt(d2)
                     d2=d2-rpi(i)
                     if (d2 < 0.0d0 ) then
               write (6,"(\"Wrongly placed grid point in gen_gridx!\")")
                     stop
                     elseif(d2 <= d2_min) then
                     d2_min=d2
                     end if
 60               continue
                  if(d2_min <= rgim) then !less than rgim
                    nld=nld+1
                    if(nld > mxlgvn)then
                      nld=nld-1
                      write(6,"('Exceeds current Dipole grid limits,',                          ' use smaller rg in option file!')")
                      stop
                    endif
                    isd (nld)=0
                    xl(1,nld)=xp(1)
                    xl(2,nld)=xp(2)
                    xl(3,nld)=xp(3)
                  endif
               endif
 40   continue
      endif
      n_inner=nld
      write(6,1003 ) nld
1003  format(/" No of inner grid dipoles              : ",i10)
      do 12 ii=0,limit_outer
         do 12 jj=0,limit_outer
            do 12 kk=0,limit_outer
               icube(ii,jj,kk)=0
               ri=ii-mid_outer
               rj=jj-mid_outer
               rk=kk-mid_outer
               d2=(ri*ri+rj*rj+rk*rk)*drg2
               if(d2 <= rg2) icube(ii,jj,kk)=1 !within rg
 12            continue
      do 21 i=i0,i1
         ipro=int(mid_outer+(xw(1,i)-center1(1))*drg_i)
         jpro=int(mid_outer+(xw(2,i)-center1(2))*drg_i)
         kpro=int(mid_outer+(xw(3,i)-center1(3))*drg_i)
         if(ipro >= 0 AND ipro <= limit_outer AND       jpro >= 0 AND jpro <= limit_outer AND       kpro >= 0 AND kpro <= limit_outer) then !within grid/dimensi
            do 31 ii=0,4
               igrid=ipro+ii-2
               do 31 jj=0,4
                  jgrid=jpro+jj-2
                  do 31 kk=0,4
                     kgrid=kpro+kk-2
                     if(icube(igrid,jgrid,kgrid) EQUAL 1)then
                        xloca=(igrid-mid_outer)*drg+center1(1)
                        yloca=(jgrid-mid_outer)*drg+center1(2)
                        zloca=(kgrid-mid_outer)*drg+center1(3)
                        ri=xloca-xw(1,i)
                        rj=yloca-xw(2,i)
                        rk=zloca-xw(3,i)
                        d2=ri*ri+rj*rj+rk*rk !distance from grid point
                        if(d2 < rp2(i)) icube(igrid,jgrid,kgrid)=0
                     endif
 31          continue
         end if
 21   continue
      no=0
      ns=0
      efmin1=0.0015d0 
      efmin2=0.0020d0 
      do 41 ii=0,limit_outer
      do 41 jj=0,limit_outer
      do 41 kk=0,limit_outer
         if(icube(ii,jj,kk) >< 0) then
            no=no+1
            xp(1)=center1(1)+(ii-mid_outer)*drg
            xp(2)=center1(2)+(jj-mid_outer)*drg
            xp(3)=center1(3)+(kk-mid_outer)*drg
            do 42 index=1,n_inner
               ri=xl(1,index)-xp(1)
               rj=xl(2,index)-xp(2)
               rk=xl(3,index)-xp(3)
               d2=ri*ri+rj*rj+rk*rk
               d2=sqrt(d2)
               if(d2 < ((drg_inner+drg)/2.0 -0.001 )) goto 41 
 42         continue
            efx=0.0d0 
            efy=0.0d0 
            efz=0.0d0 
            do 43 i=1,n_reg1
            ri=xp(1)-xw(1,i)
            rj=xp(2)-xw(2,i)
            rk=xp(3)-xw(3,i)
            d2=ri*ri+rj*rj+rk*rk
            d1=sqrt(d2)
            ddd=1.7d0 /sqrt(d1+2.0d0 )
            d3=d1 * d2
            qr=ddd*q(i)/d3
            efx=efx + qr*ri
            efy=efy + qr*rj
            efz=efz + qr*rk
43          continue
            efnorm=sqrt(efx**2+efy**2+efz**2)
            if (efnorm. lt. efmin1) goto 41 
            nld=nld+1
            if(nld > mxlgvn)then
               nld=nld-1
               write(6,"('Exceeds current Dipole grid limits,',                    ' use smaller rg in option file!')")
               stop
            endif
            if(efnorm < efmin2) then
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
      write(6,1001 ) ns
1001  format(" No of surface constrained dipoles     : ",i10)
      do i=1,nld
         rz1(i)=10000.0 
         d2_min=10000.0 
         do index=1,n_reg1
            ri=xl(1,i)-xw(1,index)
            rj=xl(2,i)-xw(2,index)
            rk=xl(3,i)-xw(3,index)
            d2=ri*ri+rj*rj+rk*rk
            rz1(i)=dmin1(rz1(i),d2)
            d2=sqrt(d2)
            d2=d2-rpi(index)
            if (d2 < 0.0d0 ) then
            write (6,"(\"Grid point, n = \",i6,\" within Vdw radius\",      \" of atom: \",i5,f10.3,\" is less than rp!\" )") i,index,d2
            stop
            elseif(d2 <= d2_min) then
            d2_min=d2
            iz(i)=index
            rz_vdw(i)=d2
            end if
         enddo
      enddo
      evdwl(ientro)=0.0 
      vdw_in=0.0d0 
      vdw_out=0.0d0 
      if ( %FALSE% ) then
      if (iprint EQUAL 1) then
      write(6,"(\"drg_i,drg,vdwsl= \",3f8.2)") drg_inner, drg, vdwsl
      write(6,"(\"n_inner, nld =  \",2i8)") n_inner, nld
      write(6,"(/,\"Contribution to VdW from individual atoms\",/)")
      end if
      do index=i0,i1
            subvdw = vdw_in + vdw_out
            r3 = rp2(index)*rpi(index)
            r6 = r3 * r3
            c6 = vdwc6(iacw(index))/rp2(index)
            xdrg=(drg_inner/3.0d0 )**3
            do i=1,n_inner
               ri=xl(1,i)-xw(1,index)
               rj=xl(2,i)-xw(2,index)
               rk=xl(3,i)-xw(3,index)
               d2=ri*ri+rj*rj+rk*rk
               d3=d2*sqrt(d2)
               d6=d2*d2*d2
            vdw_in=vdw_in+vdwsl*c6*(2.0 *r6*r3/d6/d3-3.0 *r6/d6)*xdrg
            enddo
            xdrg=(drg/3.0d0 )**3
            do i=n_inner+1,nld
               ri=xl(1,i)-xw(1,index)
               rj=xl(2,i)-xw(2,index)
               rk=xl(3,i)-xw(3,index)
               d2=ri*ri+rj*rj+rk*rk
               d3=d2*sqrt(d2)
               d6=d2*d2*d2
            vdw_out=vdw_out+vdwsl*c6*(2.0 *r6*r3/d6/d3-3.0 *r6/d6)*xdrg
            enddo
            subvdw=vdw_in+vdw_out-subvdw
      if(iprint EQUAL 1) then
      write(6,"(1x,i4,2x,a8,f10.5)")index,atom(index),subvdw
      end if
      enddo
      evdwl(ientro)=vdw_in+vdw_out
      write(6,"(/,\"Total vdw: in, out, tot: \", 3f10.5)") vdw_in,      vdw_out, evdwl(ientro)
      end if
      return
      end
