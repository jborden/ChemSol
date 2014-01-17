module chemsol
  function ef_ld (xw,q,natoms,xl,ndipole,idiel) result (da)
    !    real(8), intent(in) :: xw(3,size(xw))
    !     Electric field at lgvn dipoles is calculated from point charges.
    integer,parameter :: mxlgvn = 10000 ! this is messy, need dynamic arrays in gen_gridx
    integer,parameter :: mxatm = 500    ! this is messy, need dynamic array in main.f90 when solute is read in
    real(8),dimension(3,mxatm),intent(in) :: xw
    real(8),dimension(mxatm),intent(in) :: q
    integer,intent(in) :: natoms  ! hopefully, this will be derivable from size(q)
    real(8),dimension(3,mxlgvn),intent(in) :: xl
    integer,intent(in) :: ndipole ! hopefully, this will be derivable from size(xl)
    real(8),dimension(3,mxlgvn) :: da
    integer,intent(in) :: idiel
    real(8) :: ri,rj,rk,r2,r1,r3,ddd
    integer :: i,j
    
    !    common /reg1/ xw(3,mxatm),q(mxatm)
    !    common /scrat8/ xl(3,mxlgvn),da(3,mxlgvn)

  ! all we need, in theory, is xl, xw, q, if dynamic array allocation is used
  da = 0.d0
  do j=1,ndipole
     do  i=1,natoms
        ! get the length vector between the solute atom and point
        ri = xl(1,j)-xw(1,i)
        rj = xl(2,j)-xw(2,i)
        rk = xl(3,j)-xw(3,i)
        r2 = ri*ri + rj*rj + rk*rk
        r1 = dsqrt(r2)
        r3 = r2 * r1
        if (idiel == 1) then
           ! eq 3 in j. phys.chem.b 1997,101,5585? reported as sqrt(r1+2.0)/1.7d0 !
           ddd = 1.7d0/sqrt(r1+2.0)
           qr = ddd * q(i) / r3
        else if (idiel == 0) then
           qr = q(i) / r3
        end if
        da(1,j) = da(1,j) + qr*ri
        da(2,j) = da(2,j) + qr*rj
        da(3,j) = da(3,j) + qr*rk
     end do
  end do
  return
end function ef_ld
