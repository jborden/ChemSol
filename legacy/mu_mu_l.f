c--------------------------------------------------------------------

      subroutine mu_mu_l (n,step_,tds)
C --  Calculation of the langevin dipoles.

      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000)
      parameter (mxatm = 500)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82),
     $      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)
     $             ,q_gas(mxatm),q_mp2(mxatm)
      common /pcefaefb/efa(3,mxlgvn)
      common efal(3,mxlgvn)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),
     :                da(3,mxlgvn)
      common /pcgrid/ drg,rg,dxp0(3),
     :              rg_inner,drg_inner
      common /outer_surf/ isd (mxlgvn)
      integer*2 isd
c..................................................................
c -- Dipoles are not updated for outer surface grid points (isd(i)=1,
c    surface constrained dipoles). These are initiated in 
c    the subroutine lgvnx.

      tds = 0.0d0
      tds_sum=0.0d0
      do 10 i=1,n
          if(isd(i).eq.0) then
          efa(1,i)=efa(1,i)+efal(1,i)
          efa(2,i)=efa(2,i)+efal(2,i)
          efa(3,i)=efa(3,i)+efal(3,i)
          efna=efa(1,i)*efa(1,i)+efa(2,i)*efa(2,i)+efa(3,i)*efa(3,i)
          efna=dsqrt(efna)
          gri_sp=drg_inner
          if(i.gt.n_inner) gri_sp=drg
          call vlgvn (efna,xjunka,tds,fma,gri_sp)
          tds_sum = tds_sum + tds
          xmua(1,i)=xmua(1,i)+step_*(fma*efa(1,i)/efna-xmua(1,i))
          xmua(2,i)=xmua(2,i)+step_*(fma*efa(2,i)/efna-xmua(2,i))
          xmua(3,i)=xmua(3,i)+step_*(fma*efa(3,i)/efna-xmua(3,i))
          end if
10    continue
c     Entropy contribution:
      tds=tds_sum
c     write(6,'("ENTROPY = ",f10.3)') tds
      return
      end
