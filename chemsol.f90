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
  integer,parameter :: mxpair2 = 10000000 ! maximum amount of interacting pairs allowed ???
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
  real(8) function entropy(mass)
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
  subroutine gen_gridx (center1,nld,ientro,iflag,nvol,isd,xl,n_inner,n_reg1,rg,drg,drg_inner,rg_inner,rgim,rpi,xw,iacw, &
       vdwc6,vdwsl,atom,evdwl,rz1,iz,rz_vdw,q)
    !implicit Real*8 (a-h,o-z)
    integer,intent(inout) :: nld,nvol(mxcenter),n_inner,iz(mxlgvn)
    integer(2),intent(inout) :: isd(mxlgvn)
    integer,intent(in) :: iflag,n_reg1,ientro,iacw(mxatm)
    real(8),intent(inout) :: evdwl(mxcenter),rz1(mxlgvn),rz_vdw(mxlgvn),xl(3,mxlgvn)
    real(8),intent(in) :: rg,drg,drg_inner,rg_inner,rgim,rpi(*),xw(3,*),center1(3),vdwc6(82),vdwsl,q(*)

    character(8),intent(in) :: atom(mxatm)
    integer :: ns,i0,i1,i,limit_inner,limit_outer,mid_inner,mid_outer,no,index
    integer :: ii,jj,kk
    integer :: icube(0:133,0:133,0:133)
    integer :: igrid,jgrid,kgrid
    integer :: ipro,jpro,kpro
    real(8) :: rg2,drg2,drgi2,rgi2,rgim2,drg_inner_i,drg_i,rp2(mxatm)
    real(8) :: ri,rj,rk,d1,d2,c6,d6,r3,r6
    real(8) :: xloca,yloca,zloca
    real(8) :: xp(3)
    real(8) :: d2_min,efmin1,efmin2,efx,efy,efz,ddd
    real(8) :: d3,qr,efnorm,subvdw,vdw_in,vdw_out,xdrg
    integer :: iprint = 0 ! set this to 1 for debug information
    !     generate grid point for langevin dipole
    !     center1...center of the grid (input)
    !     nld...number of grid points (output)

    !     Set up parameters. 
    !     rg is a grid radius (outer), rgi is the same for the inner grid.
    !     drg is outer grid spacing (3A), rgim is a trim parameter
    !     for the inner grid.
    nld=0
    ns=0
    i0=1
    i1=n_reg1
    rg2=rg*rg
    drg2=drg*drg
    drgi2=drg_inner*drg_inner
    rgi2=rg_inner*rg_inner
    rgim2=rgim*rgim
    drg_inner_i=1.d0/drg_inner
    drg_i=1.0d0/drg
    do i=1,n_reg1
       rp2(i)=rpi(i)*rpi(i)
    end do

    limit_inner=int(rg_inner*2/drg_inner)
    limit_outer=int(rg*2/drg)
    !*
    !*     [make both limits into an odd number]
    !*     
    if(mod(limit_inner,2).eq.0) limit_inner=limit_inner+1
    if(mod(limit_outer,2).eq.0) limit_outer=limit_outer+1
    mid_inner=int(limit_inner/2.d0+0.5d0)
    mid_outer=int(limit_outer/2.d0+0.5d0)

    !......................................................
    !     build inner grid, if used.
    !......................................................

    if(rg_inner.gt.0.d0) then ! build inner grid
       do  ii=0,limit_inner
          do jj=0,limit_inner
             do  kk=0,limit_inner
                icube(ii,jj,kk)=0
                ri=ii-mid_inner
                rj=jj-mid_inner
                rk=kk-mid_inner
                d2=(ri*ri+rj*rj+rk*rk)*drgi2 !radius from center
                if(d2.le.rgi2) icube(ii,jj,kk)=1   !trim to a sphere
             enddo
          enddo
       end do

       !     Make cavity 
       !     Put solute on the grid nad find the nearest grid points
       !     Out of these, grid points closer than rp will be removed.
       !     The number of points removed (nvol) is proportional to the solute
       !     volume.

       nvol(ientro) = 0

       do i=i0,i1 !loop over solute atoms
          ipro=int(mid_inner+(xw(1,i)-center1(1))*drg_inner_i)
          jpro=int(mid_inner+(xw(2,i)-center1(2))*drg_inner_i)
          kpro=int(mid_inner+(xw(3,i)-center1(3))*drg_inner_i)
          if(ipro.ge.0.and.ipro.le.limit_inner.and. &
               jpro.ge.0.and.jpro.le.limit_inner.and. &
               kpro.ge.0.and.kpro.le.limit_inner) then !within grid/dimensi
             !           [check distances to 27 closest points around the
             !            nearest grid point including itself]
             do ii=1,9
                igrid=ipro+ii-5
                do jj=1,9
                   jgrid=jpro+jj-5
                   do kk=1,9
                      kgrid=kpro+kk-5
                      if(icube(igrid,jgrid,kgrid).eq.1)then
                         xloca=(igrid-mid_inner)*drg_inner+center1(1)
                         yloca=(jgrid-mid_inner)*drg_inner+center1(2)
                         zloca=(kgrid-mid_inner)*drg_inner+center1(3)
                         ri=xloca-xw(1,i)
                         rj=yloca-xw(2,i)
                         rk=zloca-xw(3,i)
                         d2=ri*ri+rj*rj+rk*rk !distance from grid point
                         if(d2.lt.rp2(i)) then
                            icube(igrid,jgrid,kgrid)=0
                            nvol(ientro) = nvol(ientro) + 1
                         end if
                      endif
                   end do
                end do
             end do
          endif
       end do

       !     Reshape grid sphere (outer region) to a solute envelope.
       do ii=0,limit_inner
          do jj=0,limit_inner
             do kk=0,limit_inner
                if(icube(ii,jj,kk).ne.0) then
                   xp(1)=center1(1)+(ii-mid_inner)*drg_inner
                   xp(2)=center1(2)+(jj-mid_inner)*drg_inner
                   xp(3)=center1(3)+(kk-mid_inner)*drg_inner
                   d2_min=10000.d0
                   !                 Get distance to a closest VdW boundary.
                   do i=i0,i1
                      ri=xp(1)-xw(1,i)
                      rj=xp(2)-xw(2,i)
                      rk=xp(3)-xw(3,i)
                      d2=ri*ri+rj*rj+rk*rk
                      d2=sqrt(d2)
                      d2=d2-rpi(i)
                      if (d2.lt.0.d0) then
                         write (6,'("Wrongly placed grid point in gen_gridx!")')
                         stop
                      elseif(d2.le.d2_min) then
                         d2_min=d2
                      end if
                   end do
                   if(d2_min.le.rgim) then !less than rgim
                      nld=nld+1
                      !......................................................................2
                      if(nld.gt.mxlgvn)then
                         nld=nld-1
                         write(6,'("Exceeds current Dipole grid limits use smaller rg in option file!")')
                         stop
                      endif
                      isd (nld)=0
                      xl(1,nld)=xp(1)
                      xl(2,nld)=xp(2)
                      xl(3,nld)=xp(3)
                   endif
                endif
             end do
          end do
       end do
    endif

    n_inner=nld
    write(6,1003) nld 
1003 format(/' No of inner grid dipoles              : ',i10)
    !......................................................
    !     build outer grid (3.0 spacing)
    !......................................................

    do ii=0,limit_outer
       do jj=0,limit_outer
          do kk=0,limit_outer
             icube(ii,jj,kk)=0
             ri=ii-mid_outer
             rj=jj-mid_outer
             rk=kk-mid_outer
             d2=(ri*ri+rj*rj+rk*rk)*drg2
             if(d2.le.rg2) icube(ii,jj,kk)=1 !within rg
          end do
       end do
    end do

    !     Make cavity 
    !     Put solute on the grid and find the nearest grid points.
    !     Out of these, grid points closer than rp will be removed.
    !     Also, grid points coordinates of which are less than 2A from
    !     the inner grid  points are removed.

    do i=i0,i1 
       !        [find distance to nearest solute atom]
       ipro=int(mid_outer+(xw(1,i)-center1(1))*drg_i)
       jpro=int(mid_outer+(xw(2,i)-center1(2))*drg_i)
       kpro=int(mid_outer+(xw(3,i)-center1(3))*drg_i)
       if(ipro.ge.0.and.ipro.le.limit_outer.and. &
            jpro.ge.0.and.jpro.le.limit_outer.and. &
            kpro.ge.0.and.kpro.le.limit_outer) then !within grid/dimensi

          do ii=0,4
             igrid=ipro+ii-2
             do jj=0,4
                jgrid=jpro+jj-2
                do kk=0,4
                   kgrid=kpro+kk-2
                   if(icube(igrid,jgrid,kgrid).eq.1)then
                      xloca=(igrid-mid_outer)*drg+center1(1)
                      yloca=(jgrid-mid_outer)*drg+center1(2)
                      zloca=(kgrid-mid_outer)*drg+center1(3)
                      ri=xloca-xw(1,i)
                      rj=yloca-xw(2,i)
                      rk=zloca-xw(3,i)
                      d2=ri*ri+rj*rj+rk*rk !distance from grid point
                      if(d2.lt.rp2(i)) icube(igrid,jgrid,kgrid)=0
                   endif
                end do
             end do
          end do
       end if
    end do

    no=0
    ns=0
    efmin1=0.0015d0
    efmin2=0.0020d0
    !     efmin1=0.0020d0
    !     efmin2=0.0024d0

    !  do 41 ii=0,limit_outer
    iloop: do ii=0,limit_outer
       !     do 41 jj=0,limit_outer
       jloop: do jj=0,limit_outer
          !        do 41 kk=0,limit_outer
          kloop: do kk=0,limit_outer
             if(icube(ii,jj,kk).ne.0) then
                no=no+1
                xp(1)=center1(1)+(ii-mid_outer)*drg
                xp(2)=center1(2)+(jj-mid_outer)*drg
                xp(3)=center1(3)+(kk-mid_outer)*drg
                !           [loop over inner grid points]
                !           adjust spacing between the inner and outer grid
                do index=1,n_inner
                   ri=xl(1,index)-xp(1)
                   rj=xl(2,index)-xp(2)
                   rk=xl(3,index)-xp(3)
                   d2=ri*ri+rj*rj+rk*rk
                   d2=sqrt(d2)
                   !if(d2.lt.((drg_inner+drg)/2.-0.001)) goto 41
                   if(d2.lt.((drg_inner+drg)/2.-0.001)) cycle kloop
                   !              if(d2.lt.drg-0.001) goto 41
                end do

                ! -- Remove grid points in the areas with small electric field from
                !    the solute. Use distance-dependent dielectric.
                ! -- For +1 charged ion, the ef>0.002  e/A**2 threshold results in the
                !    grid radius of 14.4 A.
                ! -- Finaly, gridpoints with 0.002>ef>0.0024 will be marked in the
                !    list isd (mxlgvn). Lgvn dipoles at these grid points will be kept
                !    constant in future calculations (surface constrained lgvn dipoles).

                efx=0.d0
                efy=0.d0
                efz=0.d0

                do i=1,n_reg1
                   ri=xp(1)-xw(1,i)      
                   rj=xp(2)-xw(2,i)      
                   rk=xp(3)-xw(3,i)      
                   d2=ri*ri+rj*rj+rk*rk
                   d1=sqrt(d2)
                   ddd=1.7d0/sqrt(d1+2.d0)
                   d3=d1 * d2
                   qr=ddd*q(i)/d3
                   efx=efx + qr*ri
                   efy=efy + qr*rj
                   efz=efz + qr*rk
                end do
                efnorm=sqrt(efx**2+efy**2+efz**2)
                ! if (efnorm.lt.efmin1) goto 41
                if (efnorm.lt.efmin1) cycle kloop
                nld=nld+1
                if(nld.gt.mxlgvn) then
                   nld=nld-1
                   write(6,'("Exceeds current Dipole grid limits use smaller rg in option file!")')
                   stop
                endif

                if(efnorm.lt.efmin2) then
                   isd (nld)=1
                   ns = ns + 1
                else
                   isd (nld)=0
                end if

                xl(1,nld)=xp(1)
                xl(2,nld)=xp(2)
                xl(3,nld)=xp(3)
             endif
             !41         continue
          end do kloop
       end do jloop
    end do iloop
    write(6,1001) ns
1001 format(' No of surface constrained dipoles     : ',i10) 

    !     Change the positions of outer-grid dipoles randomly in order to 
    !     eliminate artifacts related to the regularity of the cubic 
    !     grid.
    !
    !     do i = n_inner+1, nld
    !       xl(1,i) = xl(1,i) + ran2(idum)/2.d0 -0.5d0
    !       xl(2,i) = xl(2,i) + ran2(idum)/2.d0 -0.5d0
    !       xl(3,i) = xl(3,i) + ran2(idum)/2.d0 -0.5d0
    !     end do

    !     Collect distances of grid dipoles to solute nuclei (rz1) and 
    !     solute boundary (rz_vdw). Note that rz1 is a SQUARE of the dist.
    do i=1,nld
       rz1(i)=10000.0d0  
       d2_min=10000.0d0
       do index=1,n_reg1 
          ri=xl(1,i)-xw(1,index)
          rj=xl(2,i)-xw(2,index)
          rk=xl(3,i)-xw(3,index)
          d2=ri*ri+rj*rj+rk*rk
          rz1(i)=dmin1(rz1(i),d2)
          d2=sqrt(d2)
          d2=d2-rpi(index)
          if (d2.lt.0.d0) then
             write (6,'("Grid point, n = ",i6," within Vdw radius of atom: ",i5,f10.3," is less than rp!")') i,index,d2    
             stop
          elseif(d2.le.d2_min) then
             d2_min=d2
             iz(i)=index
             rz_vdw(i)=d2
          end if
       enddo
       !        write(6,'(i7,i6,f10.3)') i, iz(i), rz_vdw(i) 
    enddo

    ! --  VdW energy (9-6 formula)
    !     The minimum of the VdW curve is at the rp distance.
    !     London coeficients (vdwC6) are atom dependent, but they are not 
    !      hybridization-dependent. (e.g, (vdwc6(C(sp3)) = vdwc6(C(sp2)). 

    evdwl(ientro)=0.0
    vdw_in=0.d0
    vdw_out=0.d0
    return
  end subroutine gen_gridx

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
    integer,intent(in)    :: ndipole,iterld,iprint,n_reg1
    real(8),intent(in)    :: xl(3,mxlgvn),xw(3,*),xmua(3,mxlgvn),q_gas(*),q(*),q_mp2(*),slgvn,clgvn
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
       ! if needed, than the zan argument will be passed to the function
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
    integer,intent(in) :: nd_lgvn,n_inner,n_reg1
    real(8),intent(in) :: q(*),xw(3,*),rg_reg1,rgim
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
  subroutine lgvnx(center1,elgvn,ndipole,ientro,iterld,fsurfa,xd,da,xmua,atomfs,iz,evdwl,xl,drg_inner,drg,n_inner, &
       rz_vdw,rzcut,q,ephil1,ephil2,vdwc6,iacw,vdwsl,clgvn,slgvn,isd,rz1,atom,n_reg1,nvol,rg,rg_inner,rgim,rpi,xw)
    real(8),intent(inout) :: elgvn,fsurfa(mxcenter),xd(3,mxlgvn),da(3,mxlgvn),xmua(3,mxlgvn),atomfs(mxatm),evdwl(mxcenter)
    real(8),intent(inout) :: xl(3,mxlgvn),rz1(mxlgvn),rz_vdw(mxlgvn)
    integer,intent(inout) :: ndipole,iz(mxlgvn),nvol(mxcenter),n_inner
    integer,intent(in) :: ientro,iacw(mxatm),iterld,n_reg1
    integer(2),intent(inout) :: isd(mxlgvn)
    real(8),intent(in) :: xw(3,*),drg_inner,drg,rzcut,q(*),ephil1,ephil2,vdwc6(82),vdwsl,clgvn
    real(8),intent(in) :: center1(3),rg,rg_inner,rgim,rpi(*),slgvn
    character(8),intent(in) :: atom(mxatm)

    real(8) :: sres,tds,efn_max,vdwsur(mxatm),gri_sp,elgvn_i,efn,vlgvn_result(3),fma,elgvna,ddd,dddx,dddy,dddz,epot,rx,ry,rz, &
         rqd,fs
    integer :: idum,i,k
    elgvn=0.d0 !total lgvn energy
    !     ephil1 and ephil2 are defined in readopt, sres is surface for 
    !     large fields.

    sres= 0.0d0 !?
    !elgvn=0.d0                ! variable set to 0 again
    fsurfa(ientro)=0.d0        ! ientro = 1 at first vlgvn call
    evdwl(ientro)=0.d0         
    ndipole=0
    idum = 1
    tds = 0.0d0               ! initialize tds to 0, should be done in lgvnx
    ! after this is called, ndipole is set
    ! n_inner = 0 before call, 741 after call
    call gen_gridx(center1,ndipole,ientro,0,nvol,isd,xl,n_inner,n_reg1,rg,drg,drg_inner,rg_inner,rgim,rpi,xw,iacw, &
         vdwc6,vdwsl,atom,evdwl,rz1,iz,rz_vdw,q)

    ! --   Cartesian coordinates of point dipoles {Angstrom} are stored
    !      (it is unclear at present why 2 different variables (xl, xd)
    !      were used for coordinates of dipoles in the old pdld code.)
    do i=1,ndipole
       xd(1,i) = xl(1,i)
       xd(2,i) = xl(2,i)
       xd(3,i) = xl(3,i)
    end do


    !      Calculate magnitudes of langevin dipoles and noniterative solvation
    !      energy from the electric field scaled by the distance-dependent 
    !      dielectric constant.

    ! --   Get electric field at the positions of the dipoles (da(3,mxlgvn))
    !      call ef_ld(ndipole,1)

    !      inquire (iolength=irec) xw
    ! print the solute cartesian coordinates and charge
    ! open(file="atoms.txt",unit=46,status="replace")
    ! do i=1,n_reg1
    !    write(46,*) xw(1,i),xw(2,i),xw(3,i)
    ! end do
    ! close(46)
    ! open(file="atom-charges.txt",unit=46,status="replace")
    ! do i=1,n_reg1
    !    write(46,*) q(i)
    ! end do
    ! close(46)
    ! ! print the langevin dipoles cartesian coordinates
    ! open(file="dipoles.txt",unit=47,status="replace")
    ! do i=1,ndipole
    !    write(47,*) xl(1,i),xl(2,i),xl(3,i)
    ! end do
    ! close(47)
    da = ef_ld(xw,q,n_reg1,xl,ndipole,1)
    ! open(file="da.txt",unit=48,status="replace")
    ! do i=1,ndipole
    !    write(48,*) da(1,i),da(2,i),da(3,i)
    ! end do
    ! close(48)

    efn_max=-10.0d0
    !      elgvn = 0.d0
    do i = 1,n_reg1
       vdwsur(i) = 0.d0
    end do

    do i=1,ndipole

       gri_sp=drg_inner
       if(i.gt.n_inner) gri_sp=drg 
       if(i.eq.n_inner+1) elgvn_i = elgvn ! save elgvn_i for inner grid ld energies report

       efn=dsqrt(da(1,i)*da(1,i)+da(2,i)*da(2,i)+da(3,i)*da(3,i))
       if (efn.gt.efn_max) efn_max=efn ! set a new efn_max if a greater one is found
       ! call vlgvn(efn,elgvn,xjunk,fma,gri_sp)
       vlgvn_result = vlgvn_F(efn,gri_sp,slgvn)
       fma = vlgvn_result(1)
       tds = vlgvn_result(2)
       elgvn  =  elgvn + vlgvn_result(3)
       xmua(1,i)=fma*da(1,i)/efn
       xmua(2,i)=fma*da(2,i)/efn
       xmua(3,i)=fma*da(3,i)/efn


       ! --  Calculate hydrophobic surface
       if (i.lt.n_inner) then
          if (rz_vdw(i).le.rzcut) then

             ! --    Calculate elstat. potential at the grid point
             epot = 0.d0
             do k=1,n_reg1
                rx = xl(1,i) - xw(1,k) 
                ry = xl(2,i) - xw(2,k) 
                rz = xl(3,i) - xw(3,k) 
                rqd = sqrt (rx*rx+ry*ry+rz*rz)
                epot = epot + q(k)/rqd
             end do

             ! --    Hydrophobic energy - negative surfaces
             !      if (epot.lt.0.d0 .and. epot.gt.-ephil2) then
             !       fs = epot/ephil2 + 1.d0
             !       write(6,'(i5,2f10.5," > - ephil2, fs = ",f10.5)') iz(i),
             !    *  rz_vdw(i), epot, fs
             !      end if
             !     
             if (epot.lt.0.d0) epot=-epot
             ! --    Positive surfaces and the surface for vdw term 
             if (epot.le.ephil1) then
                fs=1.0d0
                !        write(6,'(i5,2f10.5," < ephil1, fs = ",f10.5)') iz(i),
                !    *   rz_vdw(i), epot, fs
             end if
             if((epot.gt.ephil1).and.(epot.le.ephil2)) then
                fs=1.d0-((epot-ephil1)/(ephil2-ephil1))*(1.0-sres)
                !        write(6,'(i5,2f10.5," < ephil2, fs = ",f10.5)') iz(i),
                !    *   rz_vdw(i), epot, fs
             elseif(epot.gt.ephil2) then
                fs=sres
                !        write(6,'(i5,2f10.5," > ephil2")') i, rz_vdw(i), epot
             endif
             !      endif
             atomfs(iz(i)) = atomfs(iz(i))+fs
             fsurfa(ientro)=fsurfa(ientro)+fs
             vdwsur(iz(i)) = vdwsur(iz(i)) + 1.d0
          endif
       endif
    enddo

    !     Use atom polarizabilities (vdwc6) to calculate vdW part
    !     of the solvation enthalpy.

    evdwl(ientro) = 0.d0
    do k=1,n_reg1
       evdwl(ientro) = evdwl(ientro) + vdwc6(iacw(k))*vdwsur(k)
    end do
    evdwl(ientro) = vdwsl*evdwl(ientro)

    elgvn = clgvn * elgvn
    elgvn_i = clgvn * elgvn_i
    !     write(6,'(/,"Maximum field = ",f10.4,/)') efn_max
    !------------------------------------------------------------
    write(6,1001) ndipole,elgvn, elgvn_i ! need to get here
    !-------------------------------------------------------------
    if (iterld .eq. 0) return

    !      Calculation of the initial configuration of Langevin dipoles. 
    !      (0-th step of the iterative calculation of dipole-dipole interactions).

    da = ef_ld(xw,q,n_reg1,xl,ndipole,0)

    elgvna=0.0
    !      do 20 i=1,ndipole
    do i=1,ndipole
       gri_sp=drg_inner
       if(i.gt.n_inner) gri_sp=drg
       efn=dsqrt(da(1,i)*da(1,i)+da(2,i)*da(2,i)+da(3,i)*da(3,i))
       ! call vlgvn(efn,elgvna,dumm,fma,gri_sp)
       vlgvn_result = vlgvn_F(efn,gri_sp,slgvn)
       fma = vlgvn_result(1)
       tds = vlgvn_result(2)
       elgvn  = elgvn + vlgvn_result(3)

       ! --  Dipole moments of point (langevin) dipoles are oriented along
       !     the field from solute for inner grid and outer surface dipoles. 
       !     They are oriented randomly for odd nonsurface outer grid 
       !     dipoles. (Note that this is implemented 
       !     for the zero iteration step only. For iterative lgvn,
       !     lgvn dipoles are calculated in subroutine mu_mu_l. Only those
       !     of lgvn dipoles that lie along outer surface are constrained 
       !     to be proportional to the solute field.) 

       if(i.le.n_inner .or. isd (i).eq.1 .or. mod(i,3).eq.0) then 
          ddd=3.0d0/(sqrt(rz1(i))+2.0d0)
          xmua(1,i)=ddd*fma*da(1,i)/efn
          xmua(2,i)=ddd*fma*da(2,i)/efn
          xmua(3,i)=ddd*fma*da(3,i)/efn
          !      go to 20
          cycle

          !  -- The remaining outer grid dipoles are placed randomly.
       else
          dddx=ran2(idum)
          dddy=ran2(idum)
          dddz=ran2(idum)
          ddd=sqrt(dddx*dddx+dddy*dddy+dddz*dddz)
          xmua(1,i)=dddx*fma/ddd
          xmua(2,i)=dddy*fma/ddd
          xmua(3,i)=dddz*fma/ddd
       end if
       !20    continue
    end do

    return
    !......................................................................
1001 format(' Total number of Langevin dipoles      : ',i10/ &
         ' Noniterative lgvn energy (total,inner): ',2f10.3/)
  end subroutine lgvnx
  subroutine sci_lgvn(energy,elgvni,tds,nd,icent,itl,xl,efa,efal,drg,drg_inner,ip,ip2,ip3,isd,jp,jp2, &
       n_inner,out_cut,rdcutl,slgvn,xd,xmua,da)
    real(8),intent(inout) :: tds,elgvni,efa(3,mxlgvn),efal(3,mxlgvn),xmua(3,mxlgvn),da(3,mxlgvn)
    integer,intent(inout) :: ip(0:mxlgvn),ip2(0:mxlgvn),ip3(0:mxlgvn)
    real(8),intent(in) :: xl(3,mxlgvn),drg,drg_inner,out_cut,rdcutl,slgvn,xd(3,mxlgvn)
    integer,intent(in) :: itl,nd,icent,n_inner
    integer(2),intent(inout) :: isd(mxlgvn),jp(mxpair),jp2(mxpair2)
    real(8) :: conv2,epom(44),scale,eita,eiti,eave,stepaw,energy
    integer :: i,l,nd2plot,ni,iconvp,iconv
    
    tds = 0.0d0
    write(6,1006) 
    !     write(6,1002) nd
    !     call pairlistw(nd,xd,jp3)
    !      call pairlistw(nd,xd)
    call pairlistw(nd,xd,rdcutl,out_cut,ip,ip2,ip3,isd,jp,jp2)
    write(6,1008)
    conv2=-332.d0*slgvn
    do i=1,10
       epom(i) = 0.d0
    end do
    ! --  Relaxation of the langevin dipoles
    !  do 40 l=1,itl
    do l=1,itl
       !     Dump dipoles into xyz file for viewing by the program xmol
       if (.false.) then
          open (50, file='cs_dipoles.xyz')
          nd2plot=0
          do i=1,nd
             if (xl(1,i).eq.0.0) nd2plot=nd2plot+1
          enddo
          write(50,2006) nd2plot
          write(50,*)
          do i=1,nd
             scale = 25.
             if (i.le.n_inner) scale = 50.
             if (xl(1,i).eq.0.0) then
                write (50,2007) xl(1,i),xl(2,i), xl(3,i), &
                     scale*xmua(1,i),scale*xmua(2,i),scale*xmua(3,i)
             endif
          enddo
       end if
2006   format(i5)
2007   format('X',6f10.3)


       eita=0.
       do i=1,n_inner
          eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
       enddo
       eiti=conv2*eita
       do i=n_inner+1,nd
          eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
       enddo
       eita=conv2*eita
       write(6,1013) l-1,eita,eiti,eita-tds
       !     Remember 0th iteration energy (noniterative lgvn energy).
       if (l.eq.1) elgvni=eita

       ! --  Calculate deviation of the last computed free energy from the average
       !     of previous ni iterations. If this is less than 0.1% langevin
       !     energy is converged. 
       ni=10
       if(l.gt.5) then
          eave=0.d0
          epom(ni+1)= eita - tds
          do i=1,ni
             epom(i) = epom(i+1)
             eave=eave+epom(i)
          end do

          eave=eave/float(ni)
          iconvp=1
          iconv=1
          do i=1,ni
             if((100.d0*abs((eave-epom(i))/eave)).gt.0.1d0) iconv=0
             iconvp = iconvp * iconv
          end do

          ! --  Alternatively, stop iterations if the minimum of the dG is reached.
          if(l.gt.16.and.epom(ni).gt.epom(ni-1).and. &
               epom(ni-1).gt.epom(ni-2))  iconvp = 1
          if(iconvp.eq.1) then 
             !       write(6,'("Elgvn converged. Average = ",f8.3)') eave
             ! #GOTO 50
             !goto 50
             exit
          end if
       end if

       !call newf_lcut(nd,l)
       efa = newf_lcut_f(nd,l,da,ip,jp,xd,xmua)
       ! -- In the local reaction field approximation, long-range dipole-dipole
       !    interactions are updated only in the first n steps; the field from
       !    distant dipoles is fixed at its last value afterwards.
       !     if(l.lt.10.or.l.eq.45) call updatelong(nd,l,jp3)
       if(l.lt.10.or.l.eq.45) then
          !call updatelong(nd,l)
          efal = updatelong_f(nd,l,ip2,jp2,xd,xmua)
       endif
       stepaw=0.2
       if(l.lt.4) stepaw=0.10
       if(l.lt.11) stepaw=0.20
       if(l.gt.28) stepaw=0.30
       if(l.gt.80) stepaw=0.50
       call mu_mu_l (nd,stepaw,tds,isd,efa,efal,drg_inner,n_inner,drg,xmua,slgvn)
       !40   continue
    end do

    if (iconvp.ne.1) then
       eita=0.d0
       do i=1,n_inner
          eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
       enddo
       eiti=conv2*eita
       do i=n_inner+1,nd
          eita=eita+(xmua(1,i)*da(1,i)+xmua(2,i)*da(2,i)+xmua(3,i)*da(3,i))
       enddo
       eita=conv2*eita
       write(6,1014) l-1,eita,eiti,eita-tds
    end if
    !50   continue*
    energy=eita

    return
    !....................................................................
1002 format(1x,'Total number of dipoles = ',i8)
1006 format(/,1x,'Iterative LD calculation for the chosen grid')
1008 format(//1x,'Iteration',9x,'Elgvn',5x,'Elgvn_inner',5x, &
         'Elgvn-TdS_immob'/)
1013 format(5x,i3,4x,f12.3,5x,f12.3,5x,f12.3)
1014 format(5x,i3,4x,f12.3,5x,f12.3,5x,f12.3//)
  end subroutine sci_lgvn
  subroutine dg_ld (iterld,iprint,evqdq,elgwa,evdw,elgvn,ephob,etds,ebw,clgvn,dxp0,ephil1,ephil2,iacw,ndxp, &
       pcenter,phobsl,q,q_gas,q_mp2,rg,rg_inner,rg_reg1,rgim,rpi, &
       rzcut,slgvn,tds0,vdwc6,vdwsl,xw,atom,n_reg1, &
       drg,drg_inner,rdcutl,out_cut,itl)
    real(8),intent(inout) :: ebw,elgvn,elgwa,ephob,etds,evdw,evqdq
    ! parameters in vdw.par, except for srp which seems to not be used in the rest of the program
    integer,intent(in) :: iterld,ndxp
    real(8),intent(in) :: dxp0(3),clgvn,slgvn,tds0
    ! rp is used to generate rpi
    real(8),intent(in) :: vdwc6(82)
    real(8),intent(in) :: vdwsl,phobsl,ephil1,ephil2,rzcut
    
    real(8),intent(in) :: pcenter(3),q(*),q_gas(*),q_mp2(*)
    real(8),intent(in) :: rg,rg_inner,rg_reg1,rgim,rpi(*)
    real(8),intent(in) :: xw(3,*)
    real(8),intent(in) :: drg,drg_inner,rdcutl,out_cut
    integer,intent(in) :: iprint,iacw(mxatm),n_reg1,itl
    character(8),intent(in) :: atom(*)
    real(8) :: esum,atomfs(mxatm),temp_elgvn(mxcenter),tdsl(mxcenter),tdsw_a,vatom_result(2),elgvn_ave_result(4), &
         center_new(3),da(3,mxlgvn),elgvni,efa(3,mxlgvn),efal(3,mxlgvn)
    real(8) :: oshift(3*mxcenter),rz1(mxlgvn),rz_vdw(mxlgvn),tds,temp_center(mxcenter,3)
    real(8) :: vbornx_result
    real(8) :: xd(3,mxlgvn),xl(3,mxlgvn),xmua(3,mxlgvn),fsurfa(mxcenter),evdwl(mxcenter)
    real(8) :: ecor ! used to be an input variable, not really used in the rest of the program
    integer :: i,j,iz(mxlgvn),n_inner,ndipole
    integer :: ip(0:mxlgvn),ip2(0:mxlgvn),ip3(0:mxlgvn),nvol(mxcenter)
    integer(2) :: isd(mxlgvn),jp(mxpair),jp2(mxpair2)
    esum = 0.d0
    evqdq = 0.d0
    ecor=0.d0
    do i=1,n_reg1
       atomfs(i)=0.0d0
    end do

    !      elgvn.......noniterative lgvn energy (using distance-dependent dielectric)
    !      elgvni.....lgvn energy from the 0-th iteration
    !      elgwa......iterative lgvn energy 
    !      evqdq......solute relaxation energy calculated from PCM charges
    !      ecor.......correction for electron correlation effects         

    do i=1,ndxp
       call ran_shift(i,pcenter,center_new,temp_center,ndxp,drg,drg_inner,rg_inner,dxp0,oshift)
       !        Set grid origin for the printout of dipoles into an xmol input file.
       !        center_new(1)=0.0
       !        center_new(2)=0.0
       !        center_new(3)=0.0
       !call lgvnx(center_new,elgvn,ndipole,i,iterld,fsurfa,xd,da,xmua,atomfs,iz,evdwl,xl,drg_inner,drg,n_inner,rz_vdw,rzcut,q,ephil1,ephil2,vdwc6,iacw,vdwsl,clgvn,isd,rz1,atom,n_reg1,nvol,rg,rg_inner,rgim,rpi)
       call lgvnx(center_new,elgvn,ndipole,i,iterld,fsurfa,xd,da,xmua,atomfs,iz,evdwl,xl,drg_inner,drg,n_inner, &
            rz_vdw,rzcut,q,ephil1,ephil2,vdwc6,iacw,vdwsl,clgvn,slgvn,isd,rz1,atom,n_reg1,nvol,rg,rg_inner,rgim,rpi,xw)
       esum = esum + elgvn

       if (iterld.eq.1) then
          call sci_lgvn(elgwa,elgvni,tds,ndipole,i,itl,xl,efa,efal,drg,drg_inner,ip,ip2,ip3,isd,jp,jp2,n_inner, &
               out_cut,rdcutl,slgvn,xd,xmua,da)
          temp_elgvn(i)=elgwa
          tdsl(i) = tds
       else
          elgwa= 0.d0
          temp_elgvn(i)=elgvn
          tdsl(i) = 0.d0
       end if

       tdsw_a = phobsl * fsurfa(i)
       if (iprint.eq.1) then
          write(6,'(/," Averaged atomic contributions to E_phobic:",/)')
          do j=1,n_reg1
             write(6,112) j, atom(j), phobsl*atomfs(j)/float(i)
          end do
       end if
       if (iterld.eq.1) write(6,105) elgwa,evdwl(i),tdsw_a,-tds+tdsw_a
       if (iterld.eq.0) write(6,105) elgvn,evdwl(i),tdsw_a,-tds+tdsw_a
       write(6,*) repeat("-",72)

       !call vatom (cor,vqdq,ndipole,iterld,iprint,n_reg1,xl,xw,xmua,q,q_gas,q_mp2,slgvn,clgvn)
       vatom_result = vatom_f(ndipole,iterld,iprint,n_reg1,xl,xw,xmua,q,q_gas,q_mp2,slgvn,clgvn)
       !evqdq = evqdq + vqdq
       evqdq = evqdq + vatom_result(1)
       !ecor = ecor + cor
       ecor = ecor + vatom_result(2)
    end do

    elgvn = esum/dble(ndxp)
    evqdq = - evqdq/dble(ndxp)
    ecor  = ecor/dble(ndxp)
    !call elgvn_ave(iterld,ndxp)
    !call elgvn_ave(iterld,ndxp,temp_elgvn,evdwl,fsurfa,phobsl,tdsl,nvol,evdw,ephob,etds,tds0,evqdq,elgwa)
    elgvn_ave_result = elgvn_ave_f(iterld,ndxp,temp_elgvn,evdwl,fsurfa,phobsl,tdsl,nvol,tds0,evqdq)
    evdw  = elgvn_ave_result(1)
    ephob = elgvn_ave_result(2)
    etds  = elgvn_ave_result(3)
    if (iterld.eq.1) elgwa = elgvn_ave_result(4) ! was originally in elgvn_ave
    ! tds0  = elgvn_ave_result(4)
    ! elgwa = elgvn_ave_result(5)
    !call vbornx(ndipole,ebw,center_new)
    !call vbornx(ndipole,ebw,n_inner,center_new,rg_reg1,rgim,n_reg1,q,xw)
    vbornx_result = vbornx_f(ndipole,ebw,n_inner,center_new,rg_reg1,rgim,n_reg1,q,xw)
    ebw = vbornx_result
    return
    !.......................................................................
102 format(1x,72a1)
105 format(/1x,'Elgvn = ',f8.3,'   Evdw = ',f8.3, &
         '   -TdS_phobic = ',f8.3,'   -TdS_total = ',f8.3,/)
112 format(i3,5x,a8,f9.3) 
  end subroutine dg_ld
  subroutine solvout (iterld,iprint,do_gas,evqdq,elgwa,etds,evdw,ebw,elgvn,ephob,molname,ssname)
    integer,intent(in) :: iterld ! parameter from vdw.par
    integer,intent(in) :: iprint ! parameter in main
    character(13),intent(in) :: molname ! parameter in main
    character(4),intent(in) :: ssname ! parameter in main
    logical,intent(in) :: do_gas
    ! all inout variables from 
    real(8),intent(in) :: ebw,elgvn,elgwa,ephob,etds,evdw,evqdq ! all inout variables in dg_ld, used to also have ecor as a input argument
    real(8) :: erelax,dgsolv,dgsolvni
    !      elgvn.......noniterative lgvn energy (using distance-dependent dielectric)
    !      elgwa......iterative lgvn energy 
    write(6,*) repeat ("-",72)
    write(6,*)
    write(6,*) 'FINAL RESULTS'
    !   -- Total solvation free energy from LD calculation
    erelax = evqdq
    dgsolv=elgwa+etds +evdw+ebw
    if (iterld.eq.0) dgsolv = 0.d0
    dgsolvni=elgvn+ephob+evdw+ebw
    write(6,100)
    if (iterld.eq.0) write(6,85) 
    if (iterld.eq.1) write(6,101) elgwa
    if (iterld.eq.0) write(6,101) elgvn
    write(6,102) etds 
    write(6,103) evdw
    write(6,104) ebw
    if (do_gas) write(6,105) erelax
    !     write(6,107) ecor
    if (iterld.eq.1) write(6,106) dgsolv
    if (iterld.eq.0) write(6,106) dgsolvni
    if (iprint.eq.1) then
       write(6,110)
       write(6,115) elgvn
       write(6,117) dgsolvni
    end if

    !     Short output (do not print for noniterative LD)    

    if (iterld.eq.1) then

       !     Comment the following line for IBM and uncomment the CIBM line
       open (43,file='cs.arc',access='append')
       !IBM  open (43,file='cs.arc',position='append')

       if (do_gas) then
          write(6,199)
          write(6,200) molname,ssname,elgwa, evdw, etds , erelax, ebw, &
               dgsolv-etds,dgsolv
          write(43,201) molname,ssname,elgwa, evdw, etds , erelax, ebw, &
               dgsolv-etds,dgsolv

       else
          write(6,299)
          write(6,300) molname,ssname,elgwa, evdw, etds , ebw, &
               dgsolv-etds,dgsolv
          write(43,301) molname,ssname,elgwa, evdw, etds , ebw, &
               dgsolv-etds,dgsolv
       end if
       close(43)
    end if

199 format(//,'Molecule               lgvn    VdW   -TdS   Relax',&
         '    Born   dHsolv   dGsolv')
200 format (a13,1x,a4,f9.1,3f7.1,3f9.1/)
201 format (a13,1x,a4,f9.1,3f7.1,3f9.1)
299 format(//,'Molecule               lgvn    VdW   -TdS        ',&
         '    Born   dHsolv   dGsolv')
300 format (a13,1x,a4,f9.1,2f7.1,7x,3f9.1/)
301 format (a13,1x,a4,f9.1,2f7.1,7x,3f9.1)

    return
    !......................................................................
20  FORMAT(1X,72A1//' FINAL RESULTS')
85  format(" Noniterative LD results:"/, &
         "Warning! Noniterative LD is not parametrized in cs2.1"/, &
         "         Put iterld=1 in the vdw.par file to run iterative LD"/)
100 format(//,' *****************************'/, &
         ' Langevin Dipoles Calculation '/, &
         ' *****************************'//, &
         ' Transfer of solute at 298K: GAS,1M -> WATER,1M (kcal/mol)',/)
101 format(' Langevin energy         ',10x,f20.2)
102 format(' -TdS                    ',10x,f20.2)
103 format(' VdW energy              ',10x,f20.2)
104 format(' dBorn+dOnsager          ',10x,f20.2)
105 format(' Solute relaxation       ',10x,f20.2)
107 format(' Correlation correction  ',10x,f20.2)
106 format(' Total dG                ',10x,f20.2/)

110 format(//,' ****************************'/,&
         ' Noniterative LD calculation '/,&
         ' ****************************'//,&
         ' Transfer of a solute at 298K: GAS,1M -> WATER,1M',//)
115 format(' Noniterative lgvn energy',10x,f20.2)
117 format(' dG                      ',10x,f20.2)
  end subroutine solvout
  subroutine readopt (iterld,vdwc6,dxp0,clgvn,slgvn,tds0,rp,vdwsl,phobsl,ephil1,ephil2,rzcut, &
       rpi,pcenter,rg_reg1,rg,rg_inner,rgim,ndxp,iacw,xw,q,q_gas,n_reg1,drg,iprint)
    ! parameters in vdw.par, except for srp which seems to not be used in the rest of the program
    integer,intent(inout) :: iterld,ndxp
    real(8),intent(inout) :: dxp0(3),clgvn,slgvn,tds0
    real(8),intent(inout) :: rp(82)
    real(8),intent(inout) :: vdwc6(82)
    real(8),intent(inout) :: vdwsl,phobsl,ephil1,ephil2,rzcut
    
    real(8),intent(inout) :: rpi(*),pcenter(3),rg_reg1,rg,rg_inner,rgim
    integer,intent(inout) :: iacw(mxatm)
    
    real(8),intent(in) :: xw(3,*),q(*),q_gas(*),drg
    integer,intent(in) :: n_reg1,iprint
    
    integer :: i,j,nrp,jmin
    ! srp is not used again, even though it is in vdw.par
    real(8) :: srp,r2min,r2,dmax
    character(2) :: rpinp
    !....................................................................
    !     For the adjustment of VdW radii (rp) and London coef (vdwc6) 
    !     additional input from vdw.par file is added here:
    do i=1,82
       vdwc6(i) = 0.
    end do
    read(44,1202) iterld, ndxp, dxp0(1), clgvn, slgvn,tds0
    read(44,1200) srp, rp(6), rp(7), rp(9), rp(10)
    read(44,1200) rp(13), rp(14), rp(15), rp(22), rp(24)
    read(44,1200) rp(16), rp(25)
    read(44,1200) vdwc6(1), vdwc6(6), vdwc6(9), vdwc6(13),vdwc6(22)
    read(44,1200) vdwc6(24), vdwc6(7)
    read(44,1200) vdwsl, phobsl, ephil1,ephil2, rzcut
    phobsl=phobsl/10.d0
    vdwsl=-vdwsl/10.d0

    !      Define the remaining London coef 
    do i=1, n_reg1
       if (vdwc6(iacw(i)) .eq.0.d0 ) then
          if(vdwc6(iacw(i)-1).ne.0) then
             vdwc6(iacw(i)) = vdwc6(iacw(i)-1)
          else
             vdwc6(iacw(i)) = vdwc6(iacw(i)-2)
          end if
       end if
       if (vdwc6(iacw(i)) .eq.0.d0 ) then
          if (iacw(i).le.17) vdwc6(iacw(i))=vdwc6(6)
          if (iacw(i).gt.17) vdwc6(iacw(i))=vdwc6(22)
       end if
    end do

    !     Calculate rp(H) as a linear function of the rp of the atom to        
    !     which hydrogen is covalently bonded. Note that the coefficient of 
    !     this function (srp) differs for 1st and 2nd row atoms.
    !     In addition, for inorganic oxygen (iacw=15)
    !     a separate rp(H) is used (it is read from vdw.par file).
    ! do 60 i=1,n_reg1
    do i=1,n_reg1
       if(iacw(i).ne.1 .and.iacw(i).ne.2) then
          rpi(i) = rp(iacw(i))
       else
          r2min = 100.
          jmin = 21
          !do 50 j=1,n_reg1
          do j=1,n_reg1
             !if(j.eq.i) go to 50
             if(j.eq.i) cycle
             r2=(xw(1,i)-xw(1,j))**2+(xw(2,i)-xw(2,j))**2 &
                  +(xw(3,i)-xw(3,j))**2
             if(r2.lt.r2min) then
                r2min=r2
                jmin=j
             end if
          end do
          !50         continue
          if((iacw(jmin)).eq.15) then
             iacw(i) = 2
             rpi(i) = rp(iacw(i))
             vdwc6(iacw(i)) = vdwc6(iacw(i)-1)
          else
             if(iacw(jmin).lt.18) rpi(i) = srp * rp(iacw(jmin))
             if(iacw(jmin).ge.18) rpi(i) = (srp-0.1)*rp(iacw(jmin))
          end if
       end if
       !60      continue
    end do

    ! --  Finally, allow for the change of rp parameter of any atom
    !     without changing its atom type. A new rp is specified in
    !     the end of the input file.
    read (45,'(a2)') rpinp
    if (rpinp.eq.'rp') then
       read (45,*) nrp
       do j=1,nrp
          read (45,*) i, rpi(i)
       end do
    end if

    if (iprint.eq.1) then
       write(6,'(/,"LD parameters: ",/)')
       write(6,1202) iterld, ndxp, dxp0(1), clgvn, slgvn
       write(6,1200) vdwsl, phobsl, ephil1, ephil2, rzcut
    end if

    write(6,95)
    do i=1,n_reg1
       write (6,102) i,(xw(j,i),j=1,3),q(i),q_gas(i), &
            iacw(i), rpi(i),vdwc6(iacw(i))
    enddo
    !write(6,103) dash
    write(6,*) 
    write(6,*) repeat("-",72)
    write(6,*)

95  format(//,1x,'atom #',5x,'x',8x, &
         'y',8x,'z',3x,'Q_pcm',2x,'Q_gas',1x,'atom_type', &
         3x,'rp',4x,'VdWC6'/)
102 format (i6,3f9.3,2f6.2,i6,7x,f4.2,4x,f4.2)
103 format (/1x,72a1/)
1200 format (5f6.3)
1202 format (i2,i3,4f6.3)

    ! --  Write input file for Polaris (This file is used by x_prep
    !     to create entry into amino-acid library and to create pdb file.
    !      call getenv('XPOLOUT',pname)
    !      open (36, file=pname)
    !      write(36,'(a13)') molname
    !      do i=1,n_reg1
    !      nbond=0
    !      do jj=1,10
    !      nb(jj)=0.
    !      end do
    !      do j=1,n_reg1
    !        d1=(xw(1,i)-xw(1,j))**2+(xw(2,i)-xw(2,j))**2
    !     $      +(xw(3,i)-xw(3,j))**2
    !        d1=sqrt(d1)
    !        if (i.ne.j.and.(d1.lt.1.8).and.(iacw(i).ne.1.or.
    !     $      iacw(j).ne.1)) then
    !        nbond=nbond+1
    !        nb(nbond)=j
    !        end if
    !      end do
    !      write(36,1211) i, atom(i), molname, xw(1,i), xw(2,i),
    !     $               xw(3,i), q(i), nbond, (nb(k),k=1,nbond)
    !      end do
    !1211  format(i3,1x,A4,2x,a3,4f9.3, 10i3)
    !      close (36)

    pcenter(1)=0.
    pcenter(2)=0.
    pcenter(3)=0.

    do i=1,n_reg1
       pcenter(1)=pcenter(1)+xw(1,i)/n_reg1
       pcenter(2)=pcenter(2)+xw(2,i)/n_reg1
       pcenter(3)=pcenter(3)+xw(3,i)/n_reg1
    enddo

    ! -- Find maximal radius of the solute wrt grid center.
    dmax=0.d0
    do i=1,n_reg1
       rg_reg1=(pcenter(1)-xw(1,i))**2 + (pcenter(2)-xw(2,i))**2 + &
            (pcenter(3)-xw(3,i))**2
       if (rg_reg1.gt.dmax) dmax=rg_reg1
    end do
    rg_reg1=sqrt(dmax)+2.4d0
    ! -- New grid radii are measured wrt solute surface. This ensures
    !    that enough space is attributed to the grid for large molecules.
    !    For actual choice of grid extension, other criteria are applied:
    !    The distance from the VdW surface (inner, 1A grid) and magnitude of 
    !    the field at the grid point (outer, 3A grid) - see gen_gridx subroutine.
    rg=rg_reg1+rg
    rg_inner=rg_reg1+rgim+2.d0

    if (iprint.eq.1) then
       write(6,106) pcenter
       write(6,107) rg_reg1
    end if

    if (drg.lt.3.) then
       write(6,996) 
       stop
    elseif (ndxp.gt.mxcenter) then
       write(6,997) mxcenter
       stop
    endif

106 format(/,1x,'Center of the cubic grid is defined as centroid of', &
         ' solute atoms',//,20x,'center -',3f8.3)
107 format(/,'Radius of the centroid: ',f6.2,//)
996 FORMAT(//' PROGRAM EXIT: DRG SHOULD BE AT LEAST 3.O A'//)
997 FORMAT(//' PROGRAM EXIT: MAXIMUN # OF GRIDS IS',I3//)

    return
  end subroutine readopt
end module chemsol
