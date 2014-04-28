SUBROUTINE maxnorm(n,u1,u2,out)
USE Precision
USE Constants
IMPLICIT NONE
INTEGER :: n, i
REAL(KIND=long) :: u1(n), u2(n), out
out = zero
DO i = 1, n
   out = MAX(ABS(u1(i)-u2(i)),out)
END DO
END SUBROUTINE maxnorm
