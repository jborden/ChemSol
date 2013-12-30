c--------------------------------------------------------------------
      function ran2 (idum)
      implicit Real*8 (a-h,o-z)
c     Returns a uniform random numbers between 0.0 and 1.0.
c     Set idum to any negative value to initialize or reinitialize
c     the sequence with the seed number equal -idum.
      Parameter (m=714025, ia=1366, ic=150889, rm=1.0/m)
      Integer ir(97)
      save ir, iy
      Data iff /0/
      if (idum.lt.0.or.iff.eq.0) then
        iff = 1
        idum=mod(ic-idum,m)
        do 11 j=1,97
          idum=mod(ia*idum+ic,m)
          ir(j) = idum
11      continue
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
      end
