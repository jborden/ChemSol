!--------------------------------------------------------------------
!     efn      magnitude of the electric field (e/A**2), efn=abs(da)
!     ea       Energy of the langevin dipole in external field (kcal/mol)
!              This value is changed
!     tds      
!     fma      induced dipole (component in the direction
!              of the field.)
!              This value is changed
!     gri_sp   grid spacing
!function  vlgvn (efn,ea,tds,fma,gri_sp) result (terms)
function vlgvn(efn,gri_sp) result (vlgvn_result)
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
  real(8) :: clgvn = 0.90
  real(8) :: slgvn = 0.48
  real(8), intent(in) :: efn
  real(8), intent(in) :: gri_sp
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
end function vlgvn
