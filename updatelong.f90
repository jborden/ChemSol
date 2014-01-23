subroutine updatelong(nd,l)
  !     subroutine updatelong(nd,l,jp3)
  ! --  Evaluate long-range dipole-dipole interactions to
  !     update electric field at langevin dipoles.

  implicit Real*8 (a-h,o-z)
  parameter (mxlgvn=10000)
  parameter (mxpair2=10000000)
  parameter (mxpair=5000000)
  !     parameter (mxpair3=1800000)
  common efal(3,mxlgvn)
  common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),&
       da(3,mxlgvn)
  integer*2 jp, jp2
  common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn),&
       jp2(mxpair2),ip3(0:mxlgvn)

  !:::  input vars
  integer nd
  !     integer*2 jp3(mxpair3)

  !:::  local vars
  integer i,l,npair,np,index,j
  real*8 d2,r1,r2,r3
  real*8 c5,rmu3a,rmu3ar
  !....................................................................
  do j=1,3
     do i=1,nd
        efal(j,i)=0.d0
     enddo
  enddo

  npair=0
  do i=1,nd-1
     np=ip2(i)-ip2(i-1)
     if(np.gt.0) then
        do index=1,np
           npair=npair+1
           j=jp2(npair)
           r1=xd(1,i)-xd(1,j)
           r2=xd(2,i)-xd(2,j)
           r3=xd(3,i)-xd(3,j)
           d2=r1*r1+r2*r2+r3*r3
           c5=-1.d0/(d2*d2*dsqrt(d2))
           rmu3a=3.*(r1*xmua(1,j)+r2*xmua(2,j)+r3*xmua(3,j))
           rmu3ar=3.*(r1*xmua(1,i)+r2*xmua(2,i)+r3*xmua(3,i))
           efal(1,i)=efal(1,i)+c5*(xmua(1,j)*d2-rmu3a*r1)  
           efal(2,i)=efal(2,i)+c5*(xmua(2,j)*d2-rmu3a*r2)  
           efal(3,i)=efal(3,i)+c5*(xmua(3,j)*d2-rmu3a*r3)  
           efal(1,j)=efal(1,j)+c5*(xmua(1,i)*d2-rmu3ar*r1)
           efal(2,j)=efal(2,j)+c5*(xmua(2,i)*d2-rmu3ar*r2)
           efal(3,j)=efal(3,j)+c5*(xmua(3,i)*d2-rmu3ar*r3)
        enddo
     endif
  enddo

  ! --  In the current implementation, interactions of dipoles 
  !     separated by more than out_cut are neglected.

  !     npair=0
  !     do i=1,nd-1
  !        np=ip3(i)-ip3(i-1)
  !        if(np.gt.0) then
  !           do index=1,np
  !              npair=npair+1
  !              j=jp3(npair)
  !              r1=xd(1,i)-xd(1,j)
  !              r2=xd(2,i)-xd(2,j)
  !              r3=xd(3,i)-xd(3,j)
  !              d2=r1*r1+r2*r2+r3*r3
  !              c5=-1.d0/(d2*d2*dsqrt(d2))
  !              rmu3a=3.*(r1*xmua(1,j)+r2*xmua(2,j)+r3*xmua(3,j))
  !              rmu3ar=3.*(r1*xmua(1,i)+r2*xmua(2,i)+r3*xmua(3,i))
  !              efal(1,i)=efal(1,i)+c5*(xmua(1,j)*d2-rmu3a*r1)  
  !              efal(2,i)=efal(2,i)+c5*(xmua(2,j)*d2-rmu3a*r2)  
  !              efal(3,i)=efal(3,i)+c5*(xmua(3,j)*d2-rmu3a*r3)  
  !              efal(1,j)=efal(1,j)+c5*(xmua(1,i)*d2-rmu3ar*r1)
  !              efal(2,j)=efal(2,j)+c5*(xmua(2,i)*d2-rmu3ar*r2)
  !              efal(3,j)=efal(3,j)+c5*(xmua(3,i)*d2-rmu3ar*r3)
  !            enddo
  !         endif
  !     enddo

  return
end subroutine updatelong
