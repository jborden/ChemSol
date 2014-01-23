subroutine newf_lcut(nd,l)
  ! --  Evaluation of the dipole-dipole interactions. New field
  !     at the position of the langevin dipole is calculated.
  !     Only dipoles within 2.5A-rcutl (inner, and inner-outer grid
  !     interactions) distance are used.

  implicit Real*8 (a-h,o-z)
  parameter (mxlgvn=10000)
  parameter (mxpair2=10000000)
  parameter (mxpair=5000000)
  common /pcefaefb/efa(3,mxlgvn)
  common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn), &
       da(3,mxlgvn)
  integer*2 jp, jp2
  common /biglist/ ip(0:mxlgvn),jp(mxpair),ip2(0:mxlgvn), &
       jp2(mxpair2),ip3(0:mxlgvn)

  !:::  input vars
  integer nd,l

  !:::  local vars
  integer i,i1,index,i2
  real*8 d2
  real*8 c5,rmu3a,rmu3ar
  integer npair,np
  real*8 r1,r2,r3
  !......................................................................
  ! --  da, the field at the point dipole originating from
  !     the solute molecule, remains unchanged. 
  do i=1,nd
     efa(1,i)=da(1,i)
     efa(2,i)=da(2,i)
     efa(3,i)=da(3,i)
  enddo

  ! --  Evaluation of the contribution of the dipole i2 to the
  !     electric field at the position of the dipole i1.
  !     This contribution is added to the field from the solute.
  npair=0
  do i1=1,nd-1
     np=ip(i1)-ip(i1-1)
     if(np.gt.0) then
        do index=1,np
           npair=npair+1
           i2=jp(npair)
           r1=xd(1,i1)-xd(1,i2)
           r2=xd(2,i1)-xd(2,i2)
           r3=xd(3,i1)-xd(3,i2)
           d2=r1*r1+r2*r2+r3*r3
           c5=-1./(d2*d2*dsqrt(d2))
           rmu3a=3.*(r1*xmua(1,i2)+r2*xmua(2,i2)+r3*xmua(3,i2))
           rmu3ar=3.*(r1*xmua(1,i1)+r2*xmua(2,i1)+r3*xmua(3,i1))
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

end subroutine newf_lcut
