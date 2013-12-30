
c--------------------------------------------------------------------

      subroutine vatom (cor,vqdq,ndipole,iterld,iprint)
      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000)
      parameter (mxatm=500)
      common /pcgrid/ drg,rg,dxp0(3),
     :              rg_inner,drg_inner
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82),
     :            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)
     :             ,q_gas(mxatm),q_mp2(mxatm)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),
     :                da(3,mxlgvn)
      common /lra/ clgvn, slgvn
      character*1 dash(72)
      dimension vq(mxatm)
      dimension rl(3)
      data dash/72*'-'/
      if (iprint.eq.1) then
      write(6,201) dash
      write(6,202) 
      end if

c      calculates potential from Langevin dipoles
c
c                 ->  ->    3
c      Vq = 332 * mu * r / r
c
      vqq = 0.d0
      vqdq = 0.d0
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

c     As a distance dependent field is used in noniterative LD calculation
c     Vq must be scaled by 1.7d0/dsqrt(rx+2.0d0)
             if (iterld.eq.1) then
             vq(j)=vq(j)+332*sp/(r2*r1)
             else
             vq(j)=vq(j)+332*(1.7d0/dsqrt(r1+2.0D0))*sp/(r2*r1)
             end if

 10       continue
          dq_gas = q_gas(j) - q(j)
          dq_mp2 = q_mp2(j) - q(j)
          if (iprint.eq.1) write(6,300) j, zan(j), q(j), vq(j), 
     $                     vq(j)*q(j), vq(j)*dq_gas, vq(j)*dq_mp2
          vqq = vqq + vq(j)*q(j)
          vqdq = vqdq + vq(j)*dq_gas
          cor = cor + vq(j)*dq_mp2
 20   continue
      
      vqq = slgvn*vqq
      if(iterld.eq.0) vqq = clgvn*vqq
      vqdq = 0.5d0*vqdq
      if(iterld.eq.0) vqdq = 0.8*vqdq
      if (iprint.eq.1) then
      write(6,400) vqq
      write(6,500) vqdq
      write(6,600) cor
      end if
     
 201  format(1x,72a1/,' Potentials at solute atoms [kcal/mol]',/)
 202  format(2x,'Atom #',2x,'Type',2x,'Charge',6x,'Vq     ',
     *       2x,'Vq*Q',2x,'Vq*dQ_pcm',2x,'Vq*dQ_mp2')
 300  format(2x,i5,3x,f4.1,2x,f5.2,1x,4(2x,f7.1))
 400  format(/,'Langevin energy calculated from potential at solute',
     *      ' atoms:',3x,f9.2,2x,'kcal/mol')
 500  format('Solute relaxation energy calculated from PCM charges:',
     *      '       ',f9.2,2x,'kcal/mol')
 600  format('Correction for electron correlation effects (dG_cor):',
     *      '       ',f9.2,2x,'kcal/mol')
 
      end
