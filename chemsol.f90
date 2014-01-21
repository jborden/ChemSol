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
    real(8),dimension(3) :: electrostatic_field_vector
    integer,intent(in) :: idiel
    real(8) :: ddd,scalar_component
    integer :: i,j
    ! calculate the screened electric field (eq. 2 and 3 in j. phys.chem.b 1997,101,5585)
    da = 0.d0
    do j=1,ndipole
       do  i=1,natoms
          if (idiel == 1) then
             ! screened electric field
             ! eq. 3 in j. phys.chem.b 1997,101,5585 screening factor
             ! norm2 is the euclidean vector norm, available in gfortran and ifort, cuda pgi does not have it
             ddd = sqrt(2.0d0 + norm2 ( xl(1:3,j) - xw(1:3,i) ))/1.7d0
             ! calculate the scalar component of eq. 2 in j. phys.chem.b 1997,101,5585
             scalar_component  = q(i) / (ddd * (norm2 ( xl(1:3,j) - xw(1:3,i) ) ** 3))
          else if (idiel == 0) then
             ! non-screened electric field
             scalar_component = q(i) / (norm2( xl(1:3,j) - xw(1:3,i) ) ** 3)
          end if
          ! the ith electrostatic field vector for the jth dipole
          electrostatic_field_vector = scalar_component * ( xl(1:3,j) - xw(1:3,i) )
          ! add the ith electrostatic field vector to the total electrostatic field vector of the jth dipole
          da(1:3,j) = da(1:3,j) + electrostatic_field_vector(1:3)
       end do
    end do
    return
  end function ef_ld
  function vlgvn_f(efn,gri_sp,clgvn,slgvn) result (vlgvn_result)
    ! --  Calculates the size of the projection of the induced
    !     (langevin) dipole in the direction of the electric
    !     field and its energy.
    !     efn.......magnitude of the electric field (e/A**2), efn=abs(da)
    !     xdrg......volume of the grid cell. It amounts to 27 A**3 for
    !               standard 3A grid.
    !     fma.......induced dipole (e*A) (component in the direction
    !               of the field.)
    !     ea (kcal/mol)....Energy of the langevin dipole in external field 
    !     ddd.......screening factor


    !  implicit Real*8 (a-h,o-z)
    ! common /lra/ clgvn, slgvn ! clgvn and slgvn are constants from vdw.par

    real(8), intent(in) :: efn
    real(8), intent(in) :: gri_sp
    real(8), intent(in) :: clgvn 
    real(8), intent(in) :: slgvn
    real(8), dimension(3) :: vlgvn_result

    real(8), parameter :: dipmax = 0.29d0             ! local 
    real(8), parameter :: aktm = 332.d0/0.6d0         ! local 
    real(8), parameter :: dddi = 3.0                  ! local
    real(8) :: fma,tds
    real(8) :: xdrg
    real(8) :: conv
    real(8) :: dip
    real(8) :: x,x2,x4,x6,ex,exm,algvn,corrf,ea
    xdrg = (gri_sp/3.d0)**3   ! xdrg is first determined in gen_gridx, used only in vlgvn (local)
    conv=-332.d0*slgvn        ! conv used throughout, likely placeholder var
    dip=sqrt(xdrg)*dipmax     ! local 
    !     Use old Langevin formula
    !     x=dip*efn*aktm
    x=dipmax*efn*aktm/dddi    ! x used throughout, likely placeholder var, dipmax is local
    x2=x*x                    ! local
    x4=x2*x2                  ! local
    x6=x2*x4                  ! local
    ex=exp(x)                 ! ex used throughout, likely placekeeper var
    exm=1.d0/ex               ! local
    algvn=(ex+exm)/(ex-exm) - 1.d0/x ! local 
    fma=dip*algvn             ! fma argument is set
    ! terms(1) = fma ! set the first term to fma
    !  vlgvn_result(1)  = 1
    !     tds =  - (exp(algvn**2)-1.d0)
    tds = -dlog(3.14159d0/(2.d0*acos(algvn))) ! tds is likely a placeholder, but set in elgvn_ave, lgvnx, mu_mu_l, and sci_lgvn
    corrf = 5d0*atan(x/27.d0)*exp(-x6/32.d0)  ! local
    tds = tds + corrf         !
    !     tds is larger for inner grid, i.e.
    !     tds = 1*tds for inner grid and 9*tds for outer grid
    tds = tds * 9.d0*xdrg**(2.0d0/3.0d0)
    tds = tds/17.3d0
    !  vlgvn_result(2) = tds
    ! ea = conv*fma*efn        ! ea argument is set
    !  vlgvn_result(3) = ea
    vlgvn_result = [fma,tds,conv*fma*efn]
    !     write(6,*)'ea:',ea
    return
    !.....................................................................
  end function vlgvn_f
end module chemsol
