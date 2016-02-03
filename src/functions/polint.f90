      SUBROUTINE polint(xa,ya,n,x,y,dy)
      ! Polynomial interpolation. Copied form "Numerical Recipes in
      ! FORTRAN 77"
      !
      ! Bo Terp Paulsen, botp@mek.dtu.dk

      INTEGER n,NMAX
      INTEGER, PARAMETER :: long=selected_real_kind(12,99)
      REAL(KIND=long) dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10) ! Largest anticipated value of n.
      !Given arrays xa and ya, each of length n, and given a value x, this routine
      !returns a value y, and an error estimate dy. If P (x) is the polynomial of
      !degree N − 1 such that P (xai) = yai , i = 1, . . . , n, then the returned
      !value y = P (x ).
      INTEGER i,m,ns
      REAL(KIND=long) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
    
      ns=1
      dif = abs(x-xa(1))
      DO i=1,n  !Here we find the index ns of the closest table entry,
          dift=abs(x-xa(i))

          IF (dift.lt.dif) THEN
              ns=i
              dif=dift
          ENDIF

          c(i)=ya(i) !and initialize the tableau of c’s and d’s.
          d(i)=ya(i)
      ENDDO 

      y=ya(ns) !This is the initial approximation to y.
      ns=ns-1

      DO m=1,n-1 !for each column of the tableau, we loop over the current c’s and d’s and update them.
          DO i=1,n-m
              ho=xa(i)-x
              hp=xa(i+m)-x
              w=c(i+1)-d(i)
              den=ho-hp
              IF(den.eq.0.)then
                 print*, 'failure in polint'
                 stop
              end if
              !This error can occur only if two input xa’s are (to within roundoff) identical.
              den=w/den
              d(i)=hp*den ! Here the c’s and d’s are updated.
              c(i)=ho*den
          ENDDO 
          IF (2*ns.lt.n-m)THEN
          !After each column in the tableau is completed, we decide
          !which correction, c or d, we want to add to our accu-
          !mulating value of y, i.e., which path to take through
          !the tableau—forking up or down. We do this in such a
          !way as to take the most “straight line” route through the
          !tableau to its apex, updating ns accordingly to keep track
          !of where we are. This route keeps the partial approxima-
          !tions centered (insofar as possible) on the target x. The
          !last dy added is thus the error indication.
              dy=c(ns+1)
          ELSE
              dy=d(ns)
              ns=ns-1
          ENDIF
      y=y+dy
      ENDDO 
      RETURN
      END
      

      SUBROUTINE polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
      INTEGER m,n,NMAX,MMAX
      INTEGER, PARAMETER :: long=selected_real_kind(12,99)
      REAL(KIND=long) dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      PARAMETER (NMAX=20,MMAX=20)! Maximum expected values of n and m.
    !C USES polint
      !Given arrays x1a(1:m) and x2a(1:n) of independent variables, and an m by n array
      !of function values ya(1:m,1:n), tabulated at the grid points defined by x1a and
      !x2a; and given values x1 and x2 of the independent variables; this routine
      !returns an interpolated function value y, and an accuracy indication dy 
      !(based only on the interpolation in the x1 direction, however).
      INTEGER j,k
      REAL(kind=long) ymtmp(MMAX),yntmp(NMAX)
      DO j=1,m
         DO k=1,n
            yntmp(k)=ya(j,k)
         ENDDO 
         call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
      ENDDO 
         call polint(x1a,ymtmp,m,x1,y,dy)
      return
      END
