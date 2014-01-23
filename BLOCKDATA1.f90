!--------------------------------------------------------------------
      Block Data 
      implicit Real*8 (a-h, o-z)
      PARAMETER (MXATM=500)
      parameter (mxcenter=50)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
           n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm) &
           ,q_gas(mxatm),q_mp2(mxatm) 
      common /born/ rg_reg1, rgim
      common /pcgrid/ drg,rg,dxp0(3), &
           rg_inner,drg_inner
      common /pcdipcut/ rdcutl,out_cut
      common /pctimes/ ndxp,itl,itp
      common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
      common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
      common /lra/ clgvn, slgvn
      data rg,drg,drg_inner,rdcutl,out_cut/26.0,3.0,1.0,6.1,16.0/
      data dxp0/0.0,0.0,0.0/
      data rgim/2.0/
      data ndxp,itl,itp/10,399,5/
      data rzcut/1.1/
!              H   He  Li   Be   B    C    C1  C2   N    N1
      data rp/2.3,2.3,2.15,2.00,2.40,2.65,3.0,3.25,2.65,2.85, &
!             N2   X   O   O1  O2    F    Ne  Na   Mg  Al 
           3.2,3.0,2.32,2.65,2.8,2.46,2.5,2.58,1.82,1.70, &
!             Si  P    X  S    Cl  Ar   K    Ca   Sc  Ti
           3.1,3.2,3.0,3.2,3.16,2.8,3.06,2.38,1.5,2.00, &
!              V    Cr   Mn   Fe   Co   Ni  Cu1+  Zn  Ga   Ge
           1.91,1.89,1.93,1.84,1.57,1.50,1.88,1.55,2.00,2.50, &
!             As   Se    Br   Kr   Rb  Sr    Y    Zr  Nb    Mo 
             3.00,3.00,3.44,3.00,3.25,2.70,1.75,2.00,2.00,2.00, &
!             Tc    Ru   Rh  Pd   Ag    Cd   In   Sn   Sb  Te
             2.00,2.00,2.00,2.00,2.25,1.98,2.45,2.44,3.00,3.70, &
!              I    Xe   Cs   Ba   La   Hf   Ta   W    Re   Ir 
             3.80,3.55,3.58,2.92,2.30,3.00,3.00,3.00,3.00,3.00, &
!              Pt   Au   Hg   Tl   Pb   Bi   Po   At   Rn   Fr
             3.00,1.77,1.94,2.00,2.55,3.00,3.00,3.00,3.00,3.00, &
!              Ra   Ac
             3.50,3.00/
   end Block Data
