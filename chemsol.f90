module chemsol
  implicit none
  contains
    real function entropy(mass)
    ! molecular mass 
    real(8), intent(in) :: mass
    real(8) :: Na,kB,h,pi,T,V,e,amu,m,n,S1,S2,S
    ! Avogadro's constant (number), number of molecules in a mol
    Na = 6.023d23 ! note: should be updated to 6.022
    ! Boltzmann's constant, in m^2 kg s^-2 K^-1
    kB = 1.381d-23
    ! Planck's constant, in joule seconds
    h = 6.626d-34
    ! pi
    pi = 3.14159
    ! Temperature, in Kelvin's (should have a default input parameter)
    T = 298.15
    ! ?
    V = 1d-3
    ! Euler's number
    e = 2.71828
    ! The dalton, in kg
    amu = 1.66d-27
    ! convert the atomic mass into kilograms
    m = amu * mass
    ! molarity of gas (should have a default input parameter)
    n = 1.0
    ! calculate the total amount of molecules present
    Na = Na * n 
    ! would like to have this better commented
    S1 = sqrt((2*pi*m*kB*T)**3)/h**3
    S2 = V*sqrt(e**5)/Na
    S = S1*S2
    S = dlog(S)
    S = Na*kB*S/n
    entropy = T*S/(4.18d0*1000.d0)
    return
  end function entropy
end module chemsol
