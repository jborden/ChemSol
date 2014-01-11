module chemsol
  implicit none
  contains
    real function entropy(mass)
    ! molecular mass 
    real(8), intent(in) :: mass
    real(8) :: T,V,m,n,S1,S2,S
    real(8) :: molecular_count
    ! constants
    ! note: current accepted value is 6.022    
    real(8), parameter :: Na = 6.023d23  ! Avogadro's constant (number), number of molecules in a mol
    real(8), parameter :: kB = 1.381d-23 ! Boltzmann's constant, in m^2 kg s^-2 K^-1
    real(8), parameter :: h = 6.626d-34  ! Planck's constant, in joule seconds
    real(8), parameter :: pi = 3.14159   ! pi
    real(8), parameter :: e = 2.71828    ! Euler's number
    real(8), parameter :: amu = 1.66d-27 ! The dalton, in kg
    T = 298.15     ! Temperature, in Kelvin's
    V = 1d-3       ! Volume has square ?
    ! convert the atomic mass into kilograms
    m = amu * mass
    ! molarity of gas (should be a function argument that defaults to 1)
    n = 1.0
    ! calculate the total amount of molecules present
    molecular_count = Na * n 
    ! would like to have this better commented
    S1 = sqrt((2*pi*m*kB*T)**3)/h**3
    S2 = V*sqrt(e**5)/molecular_count
    S = S1*S2
    S = dlog(S)
    S = molecular_count*kB*S/n
    entropy = T*S/(4.18d0*1000.d0)
    return
  end function entropy
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
    real(8) :: ri,rj,rk,r2,r1,r3,qr,ddd
    integer :: i,j

    !    common /reg1/ xw(3,mxatm),q(mxatm)
    !    common /scrat8/ xl(3,mxlgvn),da(3,mxlgvn)

    ! all we need, in theory, is xl, xw, q
    da = 0.d0
    do j=1,ndipole
       do  i=1,natoms
          ! get the length vector between the solute atom and point
          ri = xl(1,j)-xw(1,i)
          rj = xl(2,j)-xw(2,i)
          rk = xl(3,j)-xw(3,i)
          r2 = ri*ri + rj*rj + rk*rk
          r1 = dsqrt(r2)
          ! calculat coloumbic interaction
          r3 = r2 * r1
          if (idiel == 1) then
             ! ddd - screening factor?
             ddd = 1.7d0/sqrt(r1+2.0)
             qr = ddd * q(i) / r3
          else if (idiel == 0) then
             qr = q(i) / r3
          end if
          ! sum up the coloumbic interactions between point charges
          da(1,j) = da(1,j) + qr*ri
          da(2,j) = da(2,j) + qr*rj
          da(3,j) = da(3,j) + qr*rk
       end do
    end do
    return
  end function ef_ld
end module chemsol
