module chemsol
  implicit none
  real(8), parameter :: kB = 1.3806488d-23 ! Boltzmann's constant, in m^2 kg s^-2 K^-1    
  real(8), parameter :: Na = 6.023d23  ! Avogadro's constant (number), number of molecules in a mol
  real(8), parameter :: h = 6.626d-34  ! Planck's constant, in joule seconds
  real(8), parameter :: pi = 3.14159d0   ! pi
  real(8), parameter :: e = 2.71828    ! Euler's number
  real(8), parameter :: amu = 1.66d-27 ! The dalton, in kg
  
  integer,parameter  :: mxlgvn = 10000 ! maximum allowed langevin dipoles,  should use dynamic arrays in gen_gridx
  integer,parameter  :: mxatm = 500 ! maximum amount of atoms allowed, should be dynamic
  integer,parameter  :: mxpair  = 5000000 ! maximum amount of pairs allowed
  integer,parameter :: mxcenter = 50 ! needed by ran_shift and elgvn_ave
  integer :: iff = 0 ! a switch that tells ran2 if it has been called or not, global state
  contains
    real(8) function ran2 (idum)
      ! implicit Real*8 (a-h,o-z)
      !     Returns a uniform random numbers between 0.0 and 1.0.
      !     Set idum to any negative value to initialize or reinitialize
      !     the sequence with the seed number equal -idum.
      integer,parameter :: m = 714025, ia = 1366, ic = 150889
      real(8),parameter :: rm = 1.0 / m
      
      ! real(8),intent(in) :: idum
      integer,intent(inout)   :: idum
      ! Parameter (m=714025, ia=1366, ic=150889, rm=1.0/m)
      ! local
      integer ir(97),iy,j
      save ir, iy ! functions should not have side effects!
      ! Data iff /0/
      if (idum.lt.0.or.iff.eq.0) then
         iff = 1
         idum=mod(ic-idum,m)
         do j=1,97
            idum=mod(ia*idum+ic,m)
            ir(j) = idum
         end do
         idum=mod(ia*idum+ic,m)
         iy=idum
      end if
      j=1+(97*iy)/m
      if(j.gt.97.or.j.lt.1) then
         print*,'Problems with the random number generator, j = ', j
         stop
      end if
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j) = idum
      return
    end function ran2
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
  subroutine ran_shift(i,center1,center2,temp_center,ndxp,drg,drg_inner,rg_inner,dxp0,oshift)
    integer,intent(in) :: i,ndxp
    real(8),intent(in) :: center1(3),drg,drg_inner,rg_inner,dxp0(3)
    ! this is a good candidate for a function as these may all be packed and unpacked
    ! into one big array
    real(8),intent(inout) :: center2(3),temp_center(mxcenter,3),oshift(3*mxcenter)
    integer :: iseed,idum,dumm,kk
    real(8) :: fact,dxp(3)
    ! common /pctimes/ ndxp,itl,itp
    ! common /pcgrid/ drg,rg,dxp0(3),&
    !      rg_inner,drg_inner
    ! common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
    ! dimension oshift(3*mxcenter)
    ! dimension dxp(3), center2(3), center1(3)
    !save oshift ! don't see why this is needed

    ! --   Initialize the random number generator and
    !      generate random origin shifts for ndxp grids.
    if (i.eq.1) then
       iseed = -931
       idum = 1
       dumm = ran2(iseed)
       do kk = 1, 3*ndxp
          oshift(kk) = ran2 (idum)
       end do
    end if

    fact=drg_inner
    if(rg_inner.eq.0.d0) fact=drg
    if(i.eq.1)fact=0.0d0
    dxp(1)=fact*(1.d0-2.d0*oshift(3*i-2))
    dxp(2)=fact*(1.d0-2.d0*oshift(3*i-1))
    dxp(3)=fact*(1.d0-2.d0*oshift(3*i))

    if(i.ne.1)then
       center2(1)=center1(1)+dxp(1)+dxp0(1)
       center2(2)=center1(2)+dxp(2)+dxp0(2)
       center2(3)=center1(3)+dxp(3)+dxp0(3)
    else
       center2(1)=center1(1)+dxp0(1)
       center2(2)=center1(2)+dxp0(2)
       center2(3)=center1(3)+dxp0(3)
    endif

    temp_center(i,1)=center2(1)
    temp_center(i,2)=center2(2)
    temp_center(i,3)=center2(3)

    !      write(6,100) center1, dxp0, dxp, center2
    write(6,100) center2

    return
    !......................................................................
100 format(/ &
                                !    s ' original grid origin  ',3f9.3/
                                !    s ' original origin shift ',3f9.3/
                                !    s ' random origin shift   ',3f9.3/
         ' Grid origin            ',3f9.3)
  end subroutine ran_shift
  function ef_ld (xw,q,natoms,xl,ndipole,idiel) result (da)
    !    real(8), intent(in) :: xw(3,size(xw))
    !     Electric field at lgvn dipoles is calculated from point charges.
    !integer,parameter :: mxlgvn = 10000 ! this is messy, need dynamic arrays in gen_gridx
    !integer,parameter :: mxatm = 500    ! this is messy, need dynamic array in main.f90 when solute is read in
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
    ! integer,parameter  :: mxlgvn=10000
    ! integer,parameter  :: mxpair=5000000
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
    !integer,parameter  :: mxlgvn = 10000
    ! integer,parameter  :: mxpair = 5000000
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
    !integer,parameter  :: mxlgvn  = 10000
    ! integer,parameter  :: mxpair  = 5000000
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
  subroutine mu_mu_l (n,step_,tds,isd,efa,efal,drg_inner,n_inner,drg,xmua,slgvn)
    !integer,parameter :: mxlgvn=10000
    integer,intent(in) :: n,n_inner
    real(8),intent(in) :: step_,efal(3,mxlgvn),slgvn,drg,drg_inner
    real(8),intent(inout) :: efa(3,mxlgvn),tds,xmua(3,mxlgvn)
    integer(2),intent(in) :: isd(mxlgvn)
    !..................................................................
    ! -- Dipoles are not updated for outer surface grid points (isd(i)=1,
    !    surface constrained dipoles). These are initiated in 
    !    the subroutine lgvnx.

    ! local
    real(8),dimension(3) :: vlgvn_result
    real(8) :: tds_sum,efna,gri_sp,fma,elgvn
    integer :: i

    tds = 0.0d0
    tds_sum=0.0d0
    do i=1,n
       if(isd(i).eq.0) then
          efa(1,i)=efa(1,i)+efal(1,i)
          efa(2,i)=efa(2,i)+efal(2,i)
          efa(3,i)=efa(3,i)+efal(3,i)
          efna=efa(1,i)*efa(1,i)+efa(2,i)*efa(2,i)+efa(3,i)*efa(3,i)
          efna=dsqrt(efna)
          gri_sp=drg_inner
          if(i.gt.n_inner) gri_sp=drg
          vlgvn_result = vlgvn_F(efna,gri_sp,slgvn)
          fma = vlgvn_result(1)
          tds = vlgvn_result(2)
          elgvn  = elgvn + vlgvn_result(3)
          
          tds_sum = tds_sum + tds
          xmua(1,i)=xmua(1,i)+step_*(fma*efa(1,i)/efna-xmua(1,i))
          xmua(2,i)=xmua(2,i)+step_*(fma*efa(2,i)/efna-xmua(2,i))
          xmua(3,i)=xmua(3,i)+step_*(fma*efa(3,i)/efna-xmua(3,i))
       end if
    end do
    !     Entropy contribution:
    tds=tds_sum
    !     write(6,'("ENTROPY = ",f10.3)') tds
    return
  end subroutine mu_mu_l
  function vatom_f (ndipole,iterld,iprint,n_reg1,xl,xw,xmua,q,q_gas,q_mp2,slgvn,clgvn) result (vatom_result)
    !implicit Real*8 (a-h,o-z)
    !parameter (mxlgvn=10000)
    !parameter (mxatm=500)
    ! common /pcgrid/ drg,rg,dxp0(3), &
    !      rg_inner,drg_inner
    ! common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
    !      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm), &
    !      q_gas(mxatm),q_mp2(mxatm)
    ! common /scrat8/ xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn), &
    !      da(3,mxlgvn)
    ! common /lra/ clgvn, slgvn
    ! character*1 dash(72)
    ! dimension vq(mxatm)
    ! dimension rl(3)
    integer,intent(in)    :: ndipole,iterld,iprint,n_reg1
    real(8),intent(in)    :: xl(3,mxlgvn),xw(3,mxatm),xmua(3,mxlgvn),q_gas(mxatm),q(mxatm),q_mp2(mxatm),slgvn,clgvn
    !real(8),intent(inout) :: cor,vqdq
    integer :: i,j,k
    real(8) :: vqq,vq(mxatm),r1,r2,sp,rl(3),dq_gas,dq_mp2,vatom_result(2)
    
    ! data dash/72*'-'/
    ! if (iprint.eq.1) then
    !    write(6,201) dash
    !    write(6,202) ! there doesn't seem to be anything here, where there should be
    ! end if

    !      calculates potential from Langevin dipoles
    !
    !                 ->  ->    3
    !      Vq = 332 * mu * r / r
    !
    vqq = 0.d0
    vatom_result(1) = 0.d0
    vatom_result(2) = 0.0
    do j=1,n_reg1
       vq(j)=0.0d0
       do i=1,ndipole
          r2=0.0d0
          sp=0.0d0
          do k=1,3
             rl(k)=-xl(k,i)+xw(k,j)
             sp=sp+rl(k)*xmua(k,i)
             r2=r2+rl(k)**2
          end do
          r1=dsqrt(r2)

          !     As a distance dependent field is used in noniterative LD calculation
          !     Vq must be scaled by 1.7d0/dsqrt(rx+2.0d0)
          if (iterld.eq.1) then
             vq(j)=vq(j)+332*sp/(r2*r1)
          else
             vq(j)=vq(j)+332*(1.7d0/dsqrt(r1+2.0D0))*sp/(r2*r1)
          end if
       end do
       dq_gas = q_gas(j) - q(j)
       dq_mp2 = q_mp2(j) - q(j)
       ! not doing this as it requires additional variables to be added
       ! if need, than the zan argument will be passed to the function
       ! if (iprint.eq.1) write(6,300) j, zan(j), q(j), vq(j), &
       !      vq(j)*q(j), vq(j)*dq_gas, vq(j)*dq_mp2
       vqq = vqq + vq(j)*q(j)
       vatom_result(1) = vatom_result(1) + vq(j)*dq_gas
       vatom_result(2) = vatom_result(2) + vq(j)*dq_mp2
    end do

    vqq = slgvn*vqq
    if(iterld.eq.0) vqq = clgvn*vqq
    vatom_result(1) = 0.5d0*vatom_result(1)
    if(iterld.eq.0) vatom_result(1) = 0.8*vatom_result(1)
    if (iprint.eq.1) then
       write(6,400) vqq
       write(6,500) vatom_result(1)
       write(6,600) vatom_result(2)
    end if

201 format(1x,72a1/,' Potentials at solute atoms [kcal/mol]',/)
202 format(2x,'Atom #',2x,'Type',2x,'Charge',6x,'Vq     ', &
         2x,'Vq*Q',2x,'Vq*dQ_pcm',2x,'Vq*dQ_mp2')
300 format(2x,i5,3x,f4.1,2x,f5.2,1x,4(2x,f7.1))
400 format(/,'Langevin energy calculated from potential at solute', &
         ' atoms:',3x,f9.2,2x,'kcal/mol')
500 format('Solute relaxation energy calculated from PCM charges:', &
         '       ',f9.2,2x,'kcal/mol')
600 format('Correction for electron correlation effects (dG_cor):', &
         '       ',f9.2,2x,'kcal/mol')

  end function vatom_f
  function elgvn_ave_f(iterld,ncenter,temp_elgvn,evdwl,fsurfa,phobsl,tdsl,nvol,tds0,evqdq) result (elgvn_ave_result)
    ! implicit Real*8 (a-h,o-z)
    ! parameter (mxcenter=50)
    ! common /pcresult/ elgwa,elgvn,ebw,erelax,evdw,evqdq,ecor
    ! common/surf/ephil1,ephil2,ephob,fsurfa(mxcenter),evdwl(mxcenter)
    ! common /pcsav_center/ temp_center(mxcenter,3),temp_elgvn(mxcenter)
    ! common /vdwtds/ vdwsl,rzcut,phobsl,tdsl(mxcenter),etds,amas,tds0
    ! common /volume/ nvol(mxcenter)

    ! character*1 dash(72)
    ! data dash/72*'-'/
    !...................................................................
    integer,intent(in) :: iterld,ncenter,nvol(mxcenter)
    real(8),intent(in) :: temp_elgvn(mxcenter),evdwl(mxcenter),fsurfa(mxcenter),phobsl
    real(8),intent(in) :: tdsl(mxcenter),evqdq,tds0
    real(8) :: elgvn_ave_result(4)
    integer :: i,j,jj,nlowest
    real(8) :: eav,eboltz,enorm,temper,elowest,wlowest,expe
    real(8) :: sum_vdwl,sum_phob,sum_tds,sum_surf,sum_vol,ave_surf,vol
    eav=0.d0
    eboltz=0.d0
    enorm=0.d0
    write(6,1000) 
    temper=300.d0 ! odd, thought this was 298.15K!

    wlowest=10000.d0
    do  jj=1,ncenter
       if(temp_elgvn(jj).lt.wlowest) then
          wlowest=temp_elgvn(jj)
          nlowest=jj
       endif
    end do

    elowest=temp_elgvn(nlowest)
    do  i=1,ncenter
       write(6,1001) i,temp_elgvn(i)
       eav=eav+temp_elgvn(i)/dble(ncenter)
       expe=dexp(-(temp_elgvn(i)-elowest)/0.002d0/temper)
       eboltz=eboltz+temp_elgvn(i)*expe
       enorm=enorm+expe
    end do

    eboltz=eboltz/enorm
    write(6,1002) elowest,eav,eboltz,temper

    sum_vdwl =0.0
    sum_phob =0.0
    sum_tds  =0.0
    sum_surf =0.0
    sum_vol = 0.0
    do j=1,ncenter
       sum_vdwl = sum_vdwl + evdwl(j)
       sum_surf = sum_surf + fsurfa(j)
       sum_phob = sum_phob+phobsl*fsurfa(j)
       sum_tds = sum_tds + tdsl(j)
       sum_vol = sum_vol + dble(nvol(j))
    enddo
    ave_surf = sum_surf/dble(ncenter)
    elgvn_ave_result(1) = sum_vdwl/dble(ncenter)
    elgvn_ave_result(2) = sum_phob/dble(ncenter)
    elgvn_ave_result(3) = sum_tds/dble(ncenter)
    vol = sum_vol/dble(ncenter)


    ! --  Add solute free-volume entropy and hydrophobic entropy:
    write(6,'(/," Contributions to -TdS (kcal/mol)")' )
    write(6,'('' Change in free-volume:  '',5x,f10.3)')  tds0
    write(6,'('' Hydrophobic          :  '',5x,f10.3)')  elgvn_ave_result(2)
    write(6,'('' Dipolar saturation   :  '',5x,f10.3)')  -elgvn_ave_result(3)

    !     write(6,'(/,''Solute volume (A**3):  '',f10.3)') vol 

    elgvn_ave_result(3) = tds0 + elgvn_ave_result(2) - elgvn_ave_result(3)

    ! there is no explicit argument for print given, thus below is commented'
    ! if (iprint.eq.1) then
    !    write(6,'(/,''Average VDW energy:            '',f10.3)') evdw
    !    write(6,'(''Average -TdS (kcal/mol):       '',f10.3)') etds
    !    write(6,'(''Average solvation entropy (eu):'',f10.3)') &
    !         -1000.d0*etds/temper
    !    write(6,'(''Average relaxation energy:     '',f10.3)') evqdq
    ! end if

    elgvn_ave_result(4) = eav
    ! if (iterld.eq.1) elgwa = eav  moved to dg_ld 
    !....................................................................
1000 format(1x,' Final LD energies for different grids'// &
         '                    Elgvn     '/)
1001 format(1x,i10,2x,f10.3)
1002 format(//' lowest energy        --->  ',f10.3/ &
         ' mean energy          --->  ',f10.3/ &
         ' boltzmann <energies> --->  ',f10.3/ &
         ' for temperature      --->  ',f10.1,' K'/)
    return
  end function elgvn_ave_f
  function vbornx_f(nd_lgvn,eb,n_inner,center1,rg_reg1,rgim,n_reg1,q,xw) result (vbornx_result)
    ! implicit Real*8 (a-h,o-z)
    ! parameter (mxlgvn=5000)
    ! parameter (mxatm=500)
    ! common /pcgrid/ drg,rg,dxp0(3), &
    !      rg_inner,drg_inner
    ! common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
    !      n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm), & 
    !      q_gas(mxatm),q_mp2(mxatm)
    ! common /born/ rg_reg1, rgim
    ! character*1 dash(72)
    ! dimension center1(3)
    ! data dash/72*'-'/
    !......................................................................
    ! write(6,201) dash ! commented out for now
    integer,intent(in) :: nd_lgvn,n_inner,n_reg1
    real(8),intent(in) :: q(mxatm),xw(3,mxatm),rg_reg1,rgim
    real(8),intent(in) :: eb,center1(3)

    real(8) :: vbornx_result
    real(8) :: chqm,dqm,dqmx,dqmy,dqmz,conv,conv_dip,v_outer,b,b3,fch,eda,eba
    integer :: i
    write(*,*)
    write(*,*) repeat("-",72)
    write(*,*)
    write(*,*) 'Estimation of the bulk energy using Born equation'
    chqm = 0.d0
    dqmx = 0.d0
    dqmy = 0.d0
    dqmz = 0.d0
    conv     =  -(1.d0-1.d0/80.d0)
    !      dGdip = -166*[(2eps-2)/(2eps+1)]*dip**2/r**3
    conv_dip =  -158.d0/161.d0
    ! -- Estimation of the radius of the cavity
    v_outer = 27.0 * float(nd_lgvn - n_inner)
    b=rg_reg1+rgim
    b=(3.*v_outer*0.25/3.14159 + b*b*b)**(1./3.)
    b3=b*b*b

    !      Onsager (dipole) approximation.   

    !      Use only molecular dipole
    !      Take only a fraction of the molec. dipole because the total dipole
    !      in the cavity = molecular dipole + sum of the Langevin dipoles,
    !      i.e the effect of the Langevin dipoles is accounted for implicitly.
    !      For charged solute, take a larger fraction of the molecular
    !      dipole as Langevin dipoles are oriented to compensate the charge.
    !      Stop using Onsager for large dipoles (to prevent wrong
    !      limit in charge separation processes for large separations.
    do i=1,n_reg1
       chqm = chqm + q(i)
       dqmx=dqmx+q(i)*(xw(1,i)-center1(1))
       dqmy=dqmy+q(i)*(xw(2,i)-center1(2))
       dqmz=dqmz+q(i)*(xw(3,i)-center1(3))
    end do
    dqm = dsqrt(dqmx**2 + dqmy**2 + dqmz**2)
    if(abs(chqm).gt.1.02) fch = 1.d0
    if(abs(chqm).gt.0.02) fch = 0.1d0
    if(abs(chqm).le.0.02) fch = 0.90d0
    if(dqm.ge.5.5d0) fch=0.d0
    eda = 166.d0*conv_dip*((fch*dqm)**2)/b3
    dqm = 4.8023*dqm
    !eb=eda
    vbornx_result = eda
    eba = 0.d0

    !      Use Born formula (charge in a cavity).
    !      For charges larger then 1, scale the Born by 0.75
    !      This is done to get a correct association curve for Na+...Na+
    if(abs(chqm).lt.1.98) eba=166.d0*conv*chqm**2/b
    if(abs(chqm).ge.1.98) eba=166.d0*conv*0.75d0*chqm**2/b
    !eb=eb+eba
    vbornx_result = vbornx_result + eba
    
    write (6,200) b,dqm,chqm,eba,eda,vbornx_result
    return
    !......................................................................
200 format(//' Born radius equals - ',f9.2,1x,'A'// &
         ' molecular dipole (Debye)        - ',f9.2/ & 
         ' molecular charge                - ',f9.2/ &
         ' Born  monopole energy           - ',f9.2/ & 
         ' scaled Onsager dipole energy    - ',f9.2// & 
         ' total bulk energy difference    - ',f9.2/)
!201 format(/1x,72a1//,' Estimation of the bulk energy using Born ', & 
!         'equation')
  end function vbornx_f
end module chemsol
