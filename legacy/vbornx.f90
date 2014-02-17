subroutine vbornx(nd_lgvn,eb,center1)
  implicit Real*8 (a-h,o-z)
  parameter (mxlgvn=5000)
  parameter (mxatm=500)
  common /pcgrid/ drg,rg,dxp0(3), &
       rg_inner,drg_inner
  common /reg1/xw(3,mxatm),zan(mxatm),q(mxatm),rp(82),vdwc6(82), &
       n_inner,n_reg1,latom(mxatm),iacw(mxatm),rpi(mxatm), & 
       q_gas(mxatm),q_mp2(mxatm)
  common /born/ rg_reg1, rgim
  character*1 dash(72)
  dimension center1(3)
  data dash/72*'-'/
  !......................................................................
  write(6,201) dash
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
  eb=eda
  eba = 0.d0

  !      Use Born formula (charge in a cavity).
  !      For charges larger then 1, scale the Born by 0.75
  !      This is done to get a correct association curve for Na+...Na+
  if(abs(chqm).lt.1.98) eba=166.d0*conv*chqm**2/b
  if(abs(chqm).ge.1.98) eba=166.d0*conv*0.75d0*chqm**2/b
  eb=eb+eba

  write (6,200) b,dqm,chqm,eba,eda,eb
  return
  !......................................................................
200 format(//' Born radius equals - ',f9.2,1x,'A'// &
       ' molecular dipole (Debye)        - ',f9.2/ & 
       ' molecular charge                - ',f9.2/ &
       ' Born  monopole energy           - ',f9.2/ & 
       ' scaled Onsager dipole energy    - ',f9.2// & 
       ' total bulk energy difference    - ',f9.2/)
201 format(/1x,72a1//,' Estimation of the bulk energy using Born ', & 
       'equation')
end subroutine vbornx
