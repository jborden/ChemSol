c---------------------------------------------------------------------


      subroutine ef_ld (ndipole,idiel)
C     Electric field at lgvn dipoles is calculated from point charges.

      implicit Real*8 (a-h,o-z)
      parameter (mxlgvn=10000)
      parameter (mxatm=500)
      common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82),
     :            n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm)
     $             ,q_gas(mxatm),q_mp2(mxatm)
      common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),
     :                da(3,mxlgvn)

      if (idiel.eq.0) then
      do 5 i=1,3
      do 5 j=1,ndipole
5     da(i,j) = 0.d0
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

      if (idiel.eq.1) then
      do 15 i=1,3
      do 15 j=1,ndipole
15    da(i,j) = 0.d0
      do 20 j=1,ndipole
      do 20 i=1,n_reg1
      ri = xl(1,j)-xw(1,i)
      rj = xl(2,j)-xw(2,i)
      rk = xl(3,j)-xw(3,i)
      r2 = ri*ri + rj*rj + rk*rk
      r1 = dsqrt(r2)
      r3 = r2 * r1
      ddd = 1.7d0/sqrt(r1+2.0)
      qr = ddd * q(i) / r3
      da(1,j) = da(1,j) + qr*ri
      da(2,j) = da(2,j) + qr*rj
      da(3,j) = da(3,j) + qr*rk
20    continue
      end if

      return
      end
