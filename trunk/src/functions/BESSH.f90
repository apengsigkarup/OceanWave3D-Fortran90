      FUNCTION BESSH (N,X)

!     This subroutine calculates the third kind modified Bessel function
!     of integer order N, for any REAL X.
      IMPLICIT NONE
      COMPLEX *16 BESSH
      REAL *8 X,BESSJ, BESSY
      INTEGER N

      BESSH = BESSJ(N,X) + (0.d0,1.d0)*BESSY(N,X)
      RETURN
      END
