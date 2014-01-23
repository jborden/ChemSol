C--------------------------------------------------------------------
      Real*8 function entropy (m)
      implicit real*8 (a-h,o-z)
      real*8 NA, kB, m, n

C     Calculate the translational entropy of the ideal gas

      NA = 6.023d23
      kB = 1.381d-23
      h = 6.626d-34
      pi = 3.14159
      T = 298.15
      V = 1d-3
      e = 2.71828
      amu = 1.66d-27
c     print*, 'enter molecular mass (amu) and molarity of gas'
c     read(*,*) m,n
      n = 1.0
      m = amu*m
      NA = NA*n

      S1 = sqrt((2*pi*m*kB*T)**3)/h**3
      S2 =  V*sqrt(e**5)/NA
      S = S1*S2
      S = dlog(S)
      S = NA*kB*S/n
c     print *, 'S [cal K-1] = ', S/4.18
      entropy = T*S/(4.18d0*1000.d0)
      return
      end
