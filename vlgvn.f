c--------------------------------------------------------------------
      subroutine vlgvn (efn,ea,tds,fma,gri_sp)
C --  Calculates the size of the projection of the induced
C     (langevin) dipole in the direction of the electric
C     field and its energy.
C     efn.......magnitude of the electric field (e/A**2), efn=abs(da)
C     xdrg......volume of the grid cell. It amounts to 27 A**3 for
C               standard 3A grid.
C     fma.......induced dipole (e*A) (component in the direction
C               of the field.)
C     ea (kcal/mol)....Energy of the langevin dipole in external field 
C     ddd.......screening factor

      implicit Real*8 (a-h,o-z)
      common /lra/ clgvn, slgvn
c......................................................................
      xdrg= (gri_sp/3.d0)**3
      conv=-332.d0*slgvn
      dipmax=0.29d0
      dip=sqrt(xdrg)*dipmax
      aktm=332.d0/0.6d0
      dddi=3.0
c        Use old Langevin formula
c        x=dip*efn*aktm
      x=dipmax*efn*aktm/dddi
      x2=x*x
      x4=x2*x2
      x6=x2*x4
      ex=exp(x)
      exm=1.d0/ex
      algvn=(ex+exm)/(ex-exm) - 1.d0/x
      fma=dip*algvn
c     tds =  - (exp(algvn**2)-1.d0)
      tds = -dlog(3.14159d0/(2.d0*acos(algvn)))
      corrf = 5d0*atan(x/27.d0)*exp(-x6/32.d0)
      tds = tds + corrf
c     tds is larger for inner grid, i.e.
c     tds = 1*tds for inner grid and 9*tds for outer grid
      tds = tds * 9.d0*xdrg**(2.0d0/3.0d0)
      tds = tds/17.3d0


      
      ea=ea+conv*fma*efn
c     write(6,*)'ea:',ea
      return
c.....................................................................
      end
