module chemsol
  implicit none
  real(8), parameter :: kB = 1.3806488d-23 ! Boltzmann's constant, in m^2 kg s^-2 K^-1    
  real(8), parameter :: Na = 6.023d23  ! Avogadro's constant (number), number of molecules in a mol
  real(8), parameter :: h = 6.626d-34  ! Planck's constant, in joule seconds
  real(8), parameter :: pi = 3.14159d0   ! pi
  real(8), parameter :: e = 2.71828    ! Euler's number
  real(8), parameter :: amu = 1.66d-27 ! The dalton, in kg
  contains
    real function entropy(mass)
    ! molecular mass 
    real(8), intent(in) :: mass
    real(8) :: T,V,m,n,S1,S2,S
    real(8) :: molecular_count
    ! constants
    ! note: current accepted value is 6.022    
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
             ! norm2 is the euclidean vector norm i.e. , available in gfortran and ifort, cuda pgi does not have it
             ddd = sqrt(2.0d0 + norm2 ( xl(1:3,j) - xw(1:3,i) ))/1.7d0
             ! calculate the scalar component of eq. 2 in j. phys.chem.b 1997,101,5585
             scalar_component  = q(i) / (ddd * (norm2 ( xl(1:3,j) - xw(1:3,i) ) ** 3))
          else if (idiel == 0) then
             ! non-screened electric field
             scalar_component = q(i) / (norm2( xl(1:3,j) - xw(1:3,i) ) ** 3)
          end if
          ! the ith electric field vector for the jth dipole
          electrostatic_field_vector = scalar_component * ( xl(1:3,j) - xw(1:3,i) )
          ! add the ith electrostatic field vector to the total electrostatic field vector of the jth dipole
          da(1:3,j) = da(1:3,j) + electrostatic_field_vector(1:3)
       end do
    end do
    return
  end function ef_ld
  function vlgvn_f(efn,gri_sp,slgvn) result (vlgvn_result)
    ! --  Calculates the size of the projection of the induced
    !     (langevin) dipole in the direction of the electric
    !     field and its energy.
    !     efn.......magnitude of the electric field (e/A**2), efn=abs(da)
    !     gri_sp....grid spacing 
    !     slgvn.....scale factor for the iterative lgvn energy

    !     xdrg......volume of the grid cell. It amounts to 27 A**3 for
    !               standard 3A grid.
    !     fma.......induced dipole (e*A) (component in the direction
    !               of the field.)
    real(8), intent(in) :: efn
    real(8), intent(in) :: gri_sp
    real(8), intent(in) :: slgvn
    real(8), dimension(3) :: vlgvn_result

    real(8), parameter :: dipmax = 0.29d0             ! local 
    real(8), parameter :: aktm = 332.d0/0.6d0         ! local 
    real(8), parameter :: dddi = 3.0                  ! screening factor, eq. 7
    real(8), parameter :: k = 0.00198720              ! boltzman constant in kcal/mol/K
    real(8), parameter :: T = 298.15                      ! temperature in Kelvins
    real(8) :: fma,tds
    real(8) :: xdrg
    real(8) :: dip
    real(8) :: x,algvn,corrf,conv
    !     Use old Langevin formula
    !     x=dip*efn*aktm
    x=dipmax*efn*aktm/dddi    ! eq.7 ? j. phys.chem.b 1999,103,10284
    !x= (dipmax*efn) / (k*T*dddi)
    algvn=cosh(x)/sinh(x) - 1.d0/x     ! eq.6 j. phys.chem.b 1999,103,10284
    xdrg = (gri_sp/3.d0)**3   ! xdrg is first determined in gen_gridx, used only in vlgvn (local)
    dip=sqrt(xdrg)*dipmax     ! local 
    fma=dip*algvn             ! fma argument is set
    ! eq.8 in j. phys.chem.b 1999,103,10284 but sign is opposite and missing R constant 3.14159d0
    tds = -dlog(pi/(2.d0*acos(algvn)))
    ! eq.9 in j. phys.chem.b 1999,103,10284, but sign is opposite and missing R constant
    corrf = 5d0*atan(x/27.d0)*exp(-(x**6)/32.d0)
    ! sum of terms in parens in eq.10 in j. phys.chem.b 1999,103,10284
    tds = tds + corrf
    !     tds is larger for inner grid, i.e.
    !     tds = 1*tds for inner grid and 9*tds for outer grid
    tds = tds * 9.d0*xdrg**(2.0d0/3.0d0)
    tds = tds/17.3d0
    conv = -332.d0*slgvn
    ! conv = 1
    vlgvn_result = [fma,tds,conv*fma*efn]
    return
  end function vlgvn_f
  subroutine pairlistw(nd,dipx,rdcutl,out_cut,ip,ip2,ip3,isd,jp,jp2)
    integer,parameter  :: mxlgvn=10000
    integer,parameter  :: mxpair=5000000
    integer,parameter  :: mxpair2=10000000
    integer,intent(in) :: nd
    integer,intent(inout) :: ip(0:mxlgvn),ip2(0:mxlgvn),ip3(0:mxlgvn)
    integer(2),intent(inout) :: jp(mxpair),jp2(mxpair2)
    integer(2),intent(in)    :: isd(mxlgvn)
    real(8),intent(in) :: dipx(3,mxlgvn),rdcutl,out_cut
    real(8)            :: rdcut2,out_cut2,r1,r2,r3,d2,rcoff,rcoff2
    integer            :: npair,npair2,npair3,ntotl,i,j,mask1
    rdcut2=rdcutl*rdcutl
    out_cut2=out_cut*out_cut

    npair=0
    ip(0)=0

    npair2=0
    ip2(0)=0

    npair3=0
    ip3(0)=0

    ntotl=0

    do i=1,nd-1
       do j=i+1,nd
          if(isd(i).eq.1.and.isd(j).eq.1) cycle
          r1=dipx(1,i)-dipx(1,j)
          r2=dipx(2,i)-dipx(2,j)
          r3=dipx(3,i)-dipx(3,j)
          d2=r1*r1+r2*r2+r3*r3
          rcoff=2.5
          !  -- Switch on/off nearest neighbors for outer grid:
          !           if(i.gt.n_inner.and.j.gt.n_inner) rcoff=drg+0.2
          !  -- Switch on/off interaction of nearest neighbors belonging to 
          !     the inner and outer grid:
          !           if(i.le.n_inner.and.j.gt.n_inner) rcoff=1.9
          rcoff2=rcoff*rcoff
          mask1=d2/rcoff2
          if(mask1.ge.1) then
             ntotl=ntotl+1
             if(d2.lt.rdcut2) then
                npair=npair+1 !possible use elsewhere, but likely local
                if (npair.gt.mxpair) then
                   write(6,1000)npair,mxpair
                   stop 
                endif
                jp(npair)=j !used elsewhere
             else if(d2.ge.rdcut2.and.d2.lt.out_cut2) then
                npair2=npair2+1 !possible use elsewhere, but likely local
                if (npair2.gt.mxpair2) then
                   write(6,1002) npair2,mxpair2
                   stop 
                endif
                jp2(npair2)=j !use elsewhere
             endif
          endif
       end do
       ip(i)=npair    !used elsewhere
       ip2(i)=npair2  ! used elsewhere
       ip3(i)=npair3  ! commented out and not used elsewhere
    end do

    write(6,1001) ntotl,npair,rdcutl,npair2,out_cut, &
         npair3,out_cut

    return
    !....................................................................
1000 format(/,1x,i8,' exceeds program mxpair dimension - ',i8)
1002 format(/,1x,i8,' exceeds program mxpair2 dimension - ',i8)
    !1003 format(/,1x,i8,' exceeds program mxpair3 dimension - ',i8,/
    !   $         1x,'please use a smaller RDCUTL'//)
1001 format(/1x,'dipole pair list'/ &
         /,1x,'       total possible pairs  - ',i10,' out of ', &
         '   2.500 A cutoff', &
         /,1x,'       pairs selected        - ',i10,' within ', &
         f8.3,' A cutoff', &
         /,1x,'       pairs selected (long) - ',i10,' within ', &
         f8.3,' A cutoff', &
         /,1x,'       pairs selected (long) - ',i10,' out of ', &
         f8.3,' A cutoff')
  end subroutine pairlistw
  function newf_lcut_f(nd,l,da,ip,jp,xd,xmua) result (efa)
    integer,parameter  :: mxlgvn = 10000
    integer,parameter  :: mxpair = 5000000
    integer,intent(in) :: nd,l,ip(0:mxlgvn)
    integer(2),intent(in) :: jp(mxpair)
    real(8),intent(in) :: da(3,mxlgvn),xd(3,mxlgvn),xmua(3,mxlgvn)

    real(8),dimension(3,mxlgvn) :: efa
    integer :: npair,i,i1,np,i2,index
    real(8) :: r1,r2,r3,d2,c5,rmu3a,rmu3ar

    !--  da, the field at the point dipole originating from
    !    the solute molecule, remains unchanged. 
    do i=1,nd
       efa(1,i)=da(1,i)
       efa(2,i)=da(2,i)
       efa(3,i)=da(3,i)
    enddo
    ! efa = da
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
  end function newf_lcut_f
  function updatelong_f(nd,l,ip2,jp2,xd,xmua) result (efal)
    integer,parameter  :: mxlgvn  = 10000
    integer,parameter  :: mxpair  = 5000000
    integer,parameter  :: mxpair2 = 10000000
    integer,intent(in) :: nd,l,ip2(0:mxlgvn)
    integer(2),intent(in) :: jp2(mxpair2)
    real(8),intent(in) :: xd(3,mxlgvn),xmua(3,mxlgvn)
    ! output var
    real(8) :: efal(3,mxlgvn)
    ! local vars
    integer :: j,i,npair,i1,np,i2,index
    real(8) :: r1,r2,r3,d2,c5,rmu3a,rmu3ar

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
    return
  end function updatelong_f
  ! function mu_mu_l (n,step_,tds) result 
  ! end function mu_mu_l
end module chemsol
